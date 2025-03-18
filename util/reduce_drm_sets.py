#!/usr/bin/env python
r'''
reduce_drm_sets.py: reduce a list of simulation-ensembles to summary CSV files

usage:
  reduce_drm_sets.py [ -E ] [ -O outfile ] [ -j N ] [ -i indexfile] ENS [...]

where:
  ENS ... is a list of simulation *directories* ("ensembles"),
and:
  -O outfile gives a template (containing exactly two occurrences of %s) for
     file outputs.  This is optional - output will go to the ENS parent
     directory, according to Sandbox conventions, if not given.
     So, -O is normally not needed, but is good for debugging.
  -j N means to use N parallel workers to process the files (default 1).
     If N = 0 or 1, no parallel workers are used: fastest for our
     small workload, and helpful for debugging.
  -i indexfile names a JSON index file (by convention, s_index.json) that 
     associates experiment names with parameter values for that experiment.
     This allows output of summary information that is labeled with the
     corresponding parameter values.
     (Giving -i is usually un-needed: if Sandbox conventions are used, this
     code automatically looks for s_index.json in the ENS parent directory)
also, importantly:
  -E means to glob-expand the given ENS into ENS/*/. In this case, s_index.json
     is sought within the given ENS, and outputs are put in ENS itself. That is,
     in this case, ENS is itself the parent directory.
     NOTE: Typically -E is used.

This program rolls up the summary information already recorded (for each ENS)
in ENS/reduce_info.csv, and places the cumulative summary in a CSV file following
the given filename template. It also tabulates some yield metrics and places in
related CSV files. Overall:
  - reduce-info.csv: one-line file with overall yield and size information
  - reduce-yield.csv: same as below but with no ensemble name information
  - reduce-yield-plus.csv: one line per ENS, with summary yield information,
    ENS name, and parameter information (if s_index.json exists). Like so:

  experiment,diam,contrast,tput,iwa,chars_earth_unique,detections_earth_all,detections_earth_unique,ensemble_size
  s_YX_D7.8_iwa44_C6.36e-11_tput0.22,7.8,6.36e-11,0.22,44,7.42,20.41,10.89,100
  s_YX_D7.4_iwa64_C4.27e-10_tput0.47,7.4,4.27e-10,0.47,64,6.85,21.02,11.51,100
  s_YX_D8.0_iwa66_C3.45e-11_tput0.13,8.0,3.45e-11,0.13,66,7.10,22.43,12.77,100
  s_YX_D9.0_iwa56_C5.86e-11_tput0.54,9.0,5.86e-11,0.54,56,11.61,33.51,18.11,100
  s_YX_D7.8_iwa50_C4.88e-10_tput0.45,7.8,4.88e-10,0.45,50,16.29,45.74,25.51,100
  s_YX_D8.0_iwa38_C5.68e-10_tput0.57,8.0,5.68e-10,0.57,38,19.88,51.44,28.32,100
(Note: the "experiment" column should have been labeled "scenario", according to
our current naming conventions.)
  - reduce-yield-all.json: one JSON object per ENS, with summary information as
    in reduce-yield-plus.csv, as well as some other information such as runtime.

This program supports nested Experiments. That is, if -E is given, and if the
parent ENS contains only Experiments (.exp suffixes -- no .fam), an omnibus 
summary is generated (ENS/reduce-yield-plus.csv) that consolidates all 
sub-Experiments (ENS/*.exp/reduce-yield-plus.csv). This allows "chunking" 
large Experiments into sub-Experiments, but still generating a single yield
summary CSV. To support this, all Experiments (ENS and ENS/*.exp) must contain
their own s_index.json (with their portion of the whole). Upstream tooling
must generate these files.

The typical use cases are as follows:
- ENS list contains N basic Ensembles (single scenario yields)
  The parent ENS should end in .exp (s_index.json expected) or .fam.
  + reduce-info.csv contains the Ensemble sizing (N, and the sum of DRM 
    counts within each Ensemble), and maximum yields.
  + reduce-yield-plus.csv contains N lines, each with average 
    yields for one Ensemble, and the Ensemble name.
    If s_index.json exists (parent is an Experiment, and ends in .exp
    by convention), this file will have extra columns giving all
    parameters for that Ensemble.
- ENS list contains N Experiments (all ending in .exp, no .fam).
  + reduce-info.csv is now summarizing summaries. It contains the number
    (N) of child Experiments, and the overall number
    of Ensembles underneath, and maximal yields (which are informative).
  + reduce-yield-plus.csv will contain yields for *all* Experiments
    below (typically many more than just N, the number of Experiments),
    together with their relevant parameter values.
- ENS list contains mixed Ensemble collections (Experiments and Families,
  ending in .exp and .fam), for a total of N collections.
  + reduce-info.csv is now summarizing summaries, but that's OK. It
    contains the number (N) of child collections, and the overall number
    of Ensembles underneath, and maximal yields. In some cases, these
    maximal yields will be for completely unrelated scenarios, and
    so not very informative, but in other cases, they will be.
  + reduce-yield-plus.csv will contain N lines, one for each collection,
    with the collection name and maximal yield. 

These cases are in fact all the same. In each case, this code is loading
ENS/reduce-info.csv for each relevant ENS, and producing a one-line
"maximal" summary, and a multi-line summary file of all these reduce-info's.

Typical usage:
  util/reduce_drm_sets.py -E sims/Yokohama_Extended.fam

'''

# turmon nov 2018, feb 2022, feb 2025

import argparse
import sys
import time
import os
import copy
import csv
import json
import warnings
import shutil
from functools import partial
from collections import OrderedDict
from collections import defaultdict
#import multiprocessing.dummy
import multiprocessing as mproc
import numpy as np

##
## Globals
##

# global verbosity mode, also usable for debugging print's
VERBOSITY = 0

# number of CPUs
N_CPU = mproc.cpu_count()

########################################
###
###  Utility Functions
###
########################################

def ensure_permissions(fn):
    r'''Ensure correct permissions on the named data file.  
    We use rw-rw-r-- = 664 (octal), to allow group-write.'''
    try:
        os.chmod(fn, 0o664)
    except OSError:
        pass # e.g., don't own the file

class WorkerMap(object):
    r'''Abstracts the multiprocessing worker-pool; switches to no workers if jobs <= 1.

    This allows you to go back to ordinary single-job processing by setting the 
    number of jobs to 1.'''
    def __init__(self, jobs):
        self.jobs = jobs
        if jobs <= 1:
            # no worker pool: just use this process
            self.pool = None
            # map function is the python map() builtin
            self.map_function = map
        else:
            # the multiprocessing pool-of-workers
            self.pool = mproc.Pool(processes=jobs)
            # the map function that comes with the above
            self.map_function = self.pool.map
    def __enter__(self):
        return self.map_function
    def __exit__(self, type, value, traceback):
        if self.pool is not None:
            self.pool.terminate()

class EnsembleRun(object):
    r'''Load and summarize an Ensemble.'''
    def __init__(self, f, sim_root):
        # allow creating a dummy object so that its properties may be queried
        self.summary = None # place-holder for later
        if f is None:
            self.name = 'dummy'
            self.Nens = 0
            self.info = {}
            return
        info = self.read_sim_summary(f, sim_root)
        # Case where 'experiment' is not in the summary (reduce-info.csv)
        # Some {FOO.exp,FOO.fam}/reduce-info.csv's do not have 'experiment'
        # defined (this was an un-noticed error). The clause below fixes this.
        if 'experiment' not in info or info['experiment'] == '':
            # print(f'No Experiment: fix with: util/{args.progname} -E {f}')
            info['experiment'] = os.path.basename(f)
        # set up object state
        self.name = f
        self.info = info
        try:
            self.Nens = info['ensemble_size']
        except KeyError:
            self.Nens = 0 # (can be empty)
    
    def convert_sim_summary(self, raw_props):
        # For each property that we want to convert: (key, converter, default_value)
        # (we don't use all of these now)
        # for 'experiment':
        # (1) we must strip the left-padding because of a 
        # hack that left-pads names with a single space to control how they
        # are ultimately displayed.  This was innocuous when names were just 
        # labels, but here they are meaningful identifiers
        # (2) If it is not present, we (later) skip the record because the Ensemble
        # is not index'able without a name - the '' signals this
        prop_map = [
            ('ensemble_size', lambda s: int(float(s)), 0), # allow '0.0'
            ('runtime', str, '2000-01-01_00:00'), # actually a date
            ('simtime', str, '2000-01-01_00:01'), # actually a date
            ('user', str, 'n/a'), 
            ('experiment', str.lstrip, ''), 
            ('detections_earth_all', float, 0.0),
            ('detections_earth_unique', float, 0.0),
            ('chars_earth_unique', float, 0.0),
            ('chars_earth_strict', float, 0.0),
            # not used in the roll-up
            ('detections_unique_mean', float, 0.0),
            ('chars_unique_mean', float, 0.0), 
            ('chars_strict_mean', float, 0.0), 
            ]
        props = {}
        for key, converter, nullval in prop_map:
            props[key] = converter(raw_props.get(key, nullval))
        return props

    def read_sim_summary(self, d, sim_root):
        r'''Summarize one simulation directory (ensemble) into a dict.'''
        # file should exist: d was screened earlier to contain this file
        info_fn = os.path.join(d, 'reduce-info.csv')
        # record root of HTML summary (directory for "index.html")
        # (this may properly be an indexer function, but it fits here)
        # (nested Experiments: cannot just take URL as basename(d))
        # info_dir is the directory containing all the summarized Ensembles
        if d.startswith(sim_root + '/'):
            info_dir = d[(len(sim_root)+1):]
        else:
            info_dir = os.path.basename(d)
        if info_dir.endswith(('.fam', '.exp')):
            # if it's a .fam/.exp, index.html is at top
            base_url = info_dir
        else:
            # if it's a plain ensemble index.html is inside html/
            base_url = os.path.join(info_dir, 'html')
        #print(f'** Baseurl: {base_url}')
        #breakpoint()
        # grab the summary data for d, if possible
        try:
            with open(info_fn) as f:
                info_items = csv.DictReader(f);
                info = self.convert_sim_summary(next(info_items)) # it is a 1-line csv
                # if clause catches (and excludes) placeholder reduce-info.csv's
                # that were created by the Dakota stub that generates scripts.
                # (they are purely artifacts and have no good information)
                if len(info) < 2:
                    # print(f'\tSkipping: {info_fn}')
                    # this will later exclude the record
                    info = dict()
                else:
                    info["url"] = base_url
        except IOError:
            # it did exist earlier, but things can happen
            info = dict()
        return info

    def extract_info(self):
        r'''Extract summary info about the Ensemble.'''
        return copy.deepcopy(self.info)

    def summarize(self, econo=True):
        r'''Find the summary of the Ensemble as a dictionary held within the object.

        If econo, delete the info and keep only the summary.'''
        # this dict holds reductions for the current sim
        summary = self.extract_info()
        # delete the base data if asked
        if econo:
            self.info = None
        # keep a reference in the object
        self.summary = summary
        # also return the summary-dictionary
        return summary


# NB: Present in outer scope as a helper for load_and_reduce() below.
def outer_load_and_reduce(f, sim_root='', verb=0):
    r'''Load a sim and summarize it into a dict.

    This must be present at the outer scope of the file so it can be loaded
    by a separate process that is created by the multiprocessing module.'''
    if verb > 0:
        print('Processing <%s> in pid #%d' % (f, os.getpid()))
    ensemble = EnsembleRun(f, sim_root)
    # dictionary for this single Ensemble
    return ensemble.summarize()

class EnsembleSummary(object):
    r'''Compute, store, and dump summary information for a set of ensembles.'''
    def __init__(self, in_files, args, lazy=True):
        '''Load an ensemble of simulations.'''
        # load index file, if given
        self.load_index(args)
        # save some useful state
        self.args = args
        self.ens_files = args.ens_files
        self.Nens = len(self.ens_files) # = 0 if no ensembles

    def load_index(self, args):
        r'''Load index mapping scenario name to parameter values.

        Sets up:
          + self.index (containing parameter values, some may be vectors)
          + self.index_csv (same values, but with vectors expanded as indexed scalars)
        The latter is used for eventual output to summary .csv files.'''
        # load index file, if given
        if args.indexfile:
            print(f'{args.progname}: Loading index file.')
            with open(args.indexfile, 'r') as fp:
                index = json.load(fp)
        else:
            print(f'{args.progname}: Not using an index file.')
            index = []
        # copy index into a dict-of-dicts, EXP_NAME -> {param1:value1, param2:value2, ...}
        self.index = dict()
        # Remove these "internal bookkeeping" fields in s_index.json
        stop_fieldnames = set(('script_name', 'run_name'))
        for s in index:
            s1 = {k:v for k,v in s.items() if k not in stop_fieldnames}
            self.index[s['run_name']] = s1
        # also create a scalarized version of the index, for CSV output
        # this is the same mapping as "index", but it expands value-lists
        # into sequences of named scalars
        index_csv = dict()
        for scenario_name, d in self.index.items():
            index_csv[scenario_name] = OrderedDict()
            for param, value in d.items():
                if isinstance(value, list):
                    # expand the vector parameter "value" element-by-element
                    for inx, v1 in enumerate(value):
                        index_csv[scenario_name]['%s%d' % (param, inx+1)] = v1
                else:
                    index_csv[scenario_name][param] = value
        # save this also in the object
        self.index_csv = index_csv

    def load_and_reduce(self):
        r'''Load sim drm/spc, reduce each sim, accumulate summaries across sims.

        Each dict gives one summary statistic over that single sim.'''

        # In general, creates a pool of workers (separate unix processes)
        # (but: if jobs <= 1, it uses ordinary python map() and does no multiprocessing)
        with WorkerMap(self.args.jobs) as map_function:
            # map the load-and-reduce function over each file
            # reductions is a list of dicts containing summaries
            # (dict is empty if no valid summary existed)
            # py3: ensure the list is materialized
            reductions = list(map_function(partial(outer_load_and_reduce,
                                                   sim_root=self.args.sim_root,
                                                   verb=self.args.verbose),
                                      self.ens_files))
        # hacky fix for no-drm case
        # note: n_valid <= len(self.ens_files)
        n_valid = len([1 for r in reductions if len(r) > 0])
        if n_valid == 0:
            reductions = outer_load_and_reduce(None)
        # must save the original data for full-results tabular output
        self.reductions = [r for r in reductions if len(r) > 0]
        # re-group the above reductions across sims
        # result is a dict containing reduced data, stored as self.summary
        self.regroup_and_accum(self.reductions)

    def regroup_and_accum(self, reductions):
        r'''Accumulate various summaries across the ensemble.

        Nothing is returned: result is placed in the object state.'''
        
        # list of attributes to accumulate
        attrs = [
            'detections_earth_all', 'detections_earth_unique', 'detections_unique_mean',
            'chars_unique_mean', 'chars_earth_unique', 'chars_earth_strict',
            'ensemble_size', 
            ]

        # 1: flatten the reductions from [ens][attribute] to [attribute]
        # accum is a dictionary of lists
        accum = {}
        for attr in attrs:
            accum[attr] = []
            for r in reductions:
                accum[attr].append(r[attr])

        # prepare to get the latest simtime -- be somewhat robust if absent
        simtimes = [r['simtime'] for r in reductions if 'simtime' in r]
        simtimes.sort()

        # 2: take sum, means, max'es (whatever) of the various attributes
        #    not all computed are actually useful, or used later
        #    some QOIs can in principle have NaNs
        summary = {}
        #np.set_printoptions(precision=3) # for interactive debugging
        for attr in attrs:
            # suppress "empty slice" warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                summary[attr + '_sum'] = np.sum(accum[attr], axis=0)
                try:
                    summary[attr + '_max'] = np.nanmax(accum[attr], axis=0)
                except ValueError:
                    summary[attr + '_max'] = np.nan # fails for empties
                summary[attr + '_mean'] = np.nanmean(accum[attr], axis=0)
                # ddof=1: we're estimating the mean separately, and want Nobs=1 => no std
                summary[attr + '_std']  = np.nanstd( accum[attr], axis=0, ddof=1)
                # number of non-NaN entries in each bin of the above averages (a vector)
                N_valid = np.sum(np.isfinite(accum[attr]), axis=0)
                summary[attr + '_nEns'] = N_valid
        # This info eventually lands in reduce-info.csv
        # Note: the info here pertains only to the top-level scenario for this run, 
        # not the child scenarios that have their own experiment names, etc., and
        # that are compiled into each line of reduce-yield-plus.csv
        extra_info = dict(
            user=os.environ['USER'],
            runtime=time.strftime("%Y-%m-%d_%H:%M"),
            simtime=simtimes[-1] if simtimes else '2000-01-01_00:00',
            # experiment: actually the "scenario" name
            # we only write the last component, not the full path (TBD)
            experiment=os.path.basename(args.expt_name_readable),
            experiment_size=len(reductions), # only already-reduced ensembles
            )
        summary.update(extra_info)
        # record this cross-ensemble summary in the object
        self.summary = summary


    def dump_results_worker(self, args, extension, first_fields, otype='csv', complain=True, extras=False):
        r'''Dump reduced data to files.

        Arguments:
          otype: string, contains "csv", "json", or both; defines output file type
          complain: boolean; do we issue warnings for out-of-date s_index.json?
          extras: boolean; set True to write all gathered scenario information
        '''

        # FIXME: JSON was added later and it shows. We should generate what-to-write,
        # and then write it to CSV, JSON, or both as a second step. (See /dev/null below.)
        # But the present implementation is completely functional.

        if 'csv' in otype:
            fn = args.outfile % (extension, 'csv')
            print('\tDumping CSV to %s' % fn)
        else:
            # wasteful, yet so expedient
           fn = '/dev/null'
        # list of basic field-names in order they should be dumped
        saved_fields = [
            'chars_earth_unique',
            'detections_earth_all',
            'detections_earth_unique',
            'ensemble_size',
            ]
        # take parameter-list, if any, from first param entry in self.index_csv
        has_index = len(self.index_csv) > 0
        if has_index:
            # any key would do - we take the first
            first_key = next(iter(self.index_csv))
            param_fields = list(self.index_csv[first_key].keys())
        else:
            param_fields = []

        # fields to dump includes everything
        all_fields = first_fields + param_fields + saved_fields
        if extras:
            # toss it all in
            extra_fields = list(set(self.reductions[0].keys()) - set(all_fields))
            all_fields.extend(extra_fields)
            saved_fields.extend(extra_fields)
        # record the final dumped data so CSV and JSON match
        d_full = []
        with open(fn, 'w') as csvfile:
            # OK to have extra fields in dict-for-row
            w = csv.DictWriter(csvfile, fieldnames=all_fields, extrasaction='ignore')
            w.writeheader()
            # dictionary mapping field -> value -- everything is a scalar here
            for r in self.reductions:
                d = {key:r[key] for key in saved_fields}
                d['experiment'] = r['experiment']
                # add in the varying parameters, if they were recorded
                if has_index:
                    scenario_name = r['experiment']
                    try:
                        scenario_params = self.index_csv[scenario_name]
                    except KeyError:
                        # happens if there are ensembles in the .exp dir, that are not
                        # in s_index.json
                        if complain:
                            print(f'{args.progname}: params for {scenario_name} not in s_index.json, is it up-to-date?')
                        # all fields (e.g., params of scenarios) will be unknowns
                        scenario_params = defaultdict(lambda: "?")
                    # add in all desired params to "d"
                    for f in param_fields:
                        d[f] = scenario_params[f]
                w.writerow(d)
                # save result for the JSON
                d_full.append(d)
        if 'csv' in otype:
            ensure_permissions(fn)
        # dump the JSON version
        if 'json' in otype:
            fn = args.outfile % (extension, 'json')
            print('\tDumping JSON to %s' % fn)
            with open(fn, 'w') as jsonfile:
                json.dump(d_full, jsonfile, indent=2)
            ensure_permissions(fn)

        
    def dump_results(self, args):
        r'''Dump multi-line summary data to files

        Exports:
          - reduce-yield.csv:
                + chars, dets, Nens
                + parameter values, if index is available
          - reduce-yield-plus.csv: above, plus
                + scenario name (recorded as "experiment")
          This means that reduce-yield.csv is not in fact very useful.
          '''

        # only complain the first time
        self.dump_results_worker(args, 'yield', [], complain=True)
        self.dump_results_worker(args, 'yield-plus', ['experiment'], complain=False)
        self.dump_results_worker(args, 'yield-all', ['experiment'],
                                 otype='json', complain=False, extras=True)

    def dump_summary(self, args):
        r'''Dump overall, one-line summary data to a file: reduce-info.csv'''

        fn = args.outfile % ('info', 'csv')
        print('\tDumping to %s' % fn)
        # list of (field-name-as-output, field-name-here),
        # in order they should be dumped
        saved_field_map = [
            ('user', 'user'),
            ('runtime', 'runtime'),
            ('simtime', 'simtime'),
            ('experiment_size', 'experiment_size'),
            ('ensemble_size', 'ensemble_size_sum'),
            ('detections_earth_all', 'detections_earth_all_max'),
            ('detections_earth_unique', 'detections_earth_unique_max'),
            ('chars_earth_unique', 'chars_earth_unique_max'),
            ('chars_earth_strict', 'chars_earth_strict_max'),
            ('experiment', 'experiment'),
            ]
        with open(fn, 'w') as csvfile:
            w = csv.DictWriter(csvfile, fieldnames=[f for f,_ in saved_field_map])
            w.writeheader()
            # dictionary mapping field -> value -- everything is a scalar here
            d = {f_out:self.summary[f_have] for f_out, f_have in saved_field_map}
            w.writerow(d)
        ensure_permissions(fn)

    def dump_readme(self, args):
        r'''Dump a README document, if possible.

        Looks in the Scripts/ directory (that corresponds to the sims/... input
        that was given in the arguments), and tries to find a README type document
        for that ensemble-set. If found, the README is copied to either README.md
        or README.html in the sims/ directory.'''

        if not args.script_root:
            print(f'\tScript directory name unavailable. Continuing.')
            return
        # map is: src: dest, where ...
        #   src = README file name in Scripts/
        #   dest = README file name in sims/
        # the goal here is to have only 2 allowed file types
        # within sims/: text (as markdown), and HTML, so that
        # the html indexer does not have to try everything
        fn_map = {
            'README.md':  'README.md',
            'README.txt': 'README.md',
            'README':     'README.md',
            'README.html': 'README.html'
            }
        # filename surgery on the given outfile template
        # (we need the "naked" directory, not the template)
        outfile_dir = os.path.dirname(args.outfile % ('dummy', 'txt'))
        # attempt to find a README in Scripts/
        for fn_src, fn_dest in fn_map.items():
            fn_full = os.path.join(args.script_root, fn_src)
            if os.path.isfile(fn_full):
                fn_out = os.path.join(outfile_dir, fn_dest)
                print(f'\tDumping {fn_src} to {fn_out}')
                shutil.copyfile(fn_full, fn_out)
                ensure_permissions(fn_out)
                break
        else:
            print(f'\tNo README present. Continuing.')

    def dump_index(self, args):
        r'''Dump the index file, if possible.

        Looks in the Scripts/ directory (that corresponds to the sims/... input
        that was given in the arguments), and tries to find s_index.json
        for that ensemble-set. If found, it is copied to the sims/ directory.'''

        if not args.indexfile:
            if args.verbose:
                print(f'\tIndex file name unavailable. Continuing.')
            return False
        # can happen if you explicitly give the index file
        # not using Sandbox => don't just put it somewhere
        if not args.indexfile.startswith('Scripts/'):
            if args.verbose:
                print(f'\tIndex file {args.indexfile} not in "Scripts/". Continuing.')
            return False
        if not os.path.isfile(args.indexfile):
            print(f'\tIndex file {args.indexfile} no longer found. Continuing.')
            return False
        #
        fn_out = args.indexfile.replace('Scripts/', 'sims/', 1)
        print(f'\tDumping {args.indexfile} to {fn_out}')
        shutil.copyfile(args.indexfile, fn_out)
        ensure_permissions(fn_out)
        return True


def main(args):
    ensemble_set = EnsembleSummary(args.infile, args)
    if ensemble_set.Nens == 0:
        print(f'{args.progname}: Warning: no actual Ensembles present.')
    else:
        print(f'{args.progname}: Found {ensemble_set.Nens} summarized ensembles.')
    # save a reference to the universe (stars and their attributes)
    print(f'{args.progname}: Reducing.')
    ensemble_set.load_and_reduce()
    print(f'{args.progname}: Dumping result tabulation.')
    ensemble_set.dump_results(args)
    print(f'{args.progname}: Dumping capsule summary.')
    ensemble_set.dump_summary(args)
    print(f'{args.progname}: Seeking README.')
    ensemble_set.dump_readme(args)
    if ensemble_set.dump_index(args):
        print(f'{args.progname}: Dumping index file.')


def expand_infile_arg(args):
    '''Expand args.infile into a list of scenarios

    Assumes args.infile has at least one member.'''
    # list of ensemble directories: our main goal
    ens_files = []
    # is this nested Experiment (foo.exp/{bar,baz,brat}.exp)?
    # (if so, we will descend into its child experiments)
    inner_experiment = False
    if args.expand:
        infiles = args.infile
        if args.infile[0].endswith('.exp'):
            # look for possible experiment-of-experiments
            # Note: not using glob.glob() because of possible efficiency issues
            # with large scenario counts - only using os.scandir() here and below
            inner_exps = []
            found_fam = False
            for entry in os.scandir(args.infile[0]):
                if entry.path.endswith('.fam'):
                    # cannot be an exp-of-exp's
                    found_fam = True
                    break
                elif entry.path.endswith('.exp'):
                    inner_exps.append(entry.path)
                else:
                    pass # don't care about other dirs
            # if Nested Experiment, replace infile X.exp with X.exp/*.exp
            if inner_exps and not found_fam:
                inner_experiment = True
                infiles = inner_exps
                print(f"{args.progname}: Found {len(infiles)} sub-experiment(s).")
        # expand into: [d/reduce-info.csv for d in infiles if is_dir(d)]
        # one-liner is possible, but want to allow granular operation timing
        for infile1 in infiles:
            for entry in os.scandir(infile1):
                if entry.is_dir():
                    if os.path.isfile(os.path.join(entry.path, 'reduce-info.csv')):
                        ens_files.append(entry.path)
    else:
        # filter: only retain input directories (members of args.infile) with reduced DRM-sets
        for f in args.infile:
            if os.path.isfile(os.path.join(f, 'reduce-info.csv')):
                ens_files.append(f)
    # one of the scenarios we're summarizing - useful later
    if not ens_files:
        args.ens1 = ''
    elif inner_experiment:
        # need the enclosing experiment
        args.ens1 = os.path.dirname(os.path.normpath(ens_files[0]))
    else:
        args.ens1 = os.path.normpath(ens_files[0])
    return ens_files


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Summarize an Experiment/Family, that is, a list of Ensembles.",
                                     epilog='Without -E, typically invoke with wildcard to look in each ENS/reduce-info.csv. With -E, typically give without wildcard, and thus look in ENS/*/reduce-info.csv')
    parser.add_argument('ens', metavar='ENS', nargs='+', default=[],
                            help='ensemble list, or (with -E) directory-of-ensembles')
    parser.add_argument('-E', '--expand', dest='expand', action='store_true', help='expand subdirectories in ENS')
    parser.add_argument('-O', '--outfile', type=str, default='',
                            help='output file template (omit to assume Sandbox conventions)')
    parser.add_argument('-i', '--indexfile', type=str, default='',
                            help='explicitly give JSON index filename (do not use with -I)')
    parser.add_argument('-I', '--index_off', dest='index_off', action='store_true', help='do not look for index file')
    parser.add_argument('-v', help='verbosity', action='count', dest='verbose',
                            default=0)
    parser.add_argument('-j', '--jobs', help=' # parallel jobs, default = %(default)d',
                      type=int, dest='jobs', default=1)
    args = parser.parse_args()
    
    # measure elapsed time
    time_start = time.time()

    # set umask in hopes that files will be group-writable
    os.umask(0o002)

    VERBOSITY = args.verbose

    # program name, for convenience
    args.progname = os.path.basename(sys.argv[0])
    
    args.infile = args.ens
    if len(args.infile) == 0:
        print('%s: Error: Need at least one input ensemble.' % args.progname)
        sys.exit(1)

    ## Expansion of file arguments
    # 1/ Unpack the input scenario-list into "args.ens_files"
    #    a/ if -E, we unpack into scenario dirs with reduce-info.csv's
    #    b/ no -E, we ensure each entry of the scenario list has reduce-info.csv
    #    c/ Generate args.ens1 = args.ens_files[0], or a dummy
    #    (this is in a dedicated function)
    # 2/ Generate other misc. names from args.ens1
    #    a/ expt_name
    #    b/ script_root
    # 3/ Establish the index filename, if any
    #    a/ guess the filename using args.ens1
    #    b/ but if -i fn, check it exists (not present for .fam's, and that's OK)
    #    c/ None is also OK (-I)
    # 4/ Get outfile
    #    a/ if explicit with -O, check it
    #    b/ if implicit, generate from expt_name
    #    c/ In either case, ensure enclosing dir exists

    # 1/
    # finalize file list (ens_files) and args.infile0
    print(f'{args.progname}: Examining {len(args.infile)} input file pattern(s) for valid scenarios.')
    args.ens_files = expand_infile_arg(args)
    print(f"{args.progname}: Found {len(args.ens_files)} ensemble summary file(s).")
    if len(args.ens_files) == 0:
        # not sure this is truly necessary, but enforcing it for now
        print(f"{args.progname}: Error: Need at least one input ensemble.")
        sys.exit(1)

    # 2/ misc filenames
    # get the experiment name from the directory - brittle, but expedient.
    # examples of operation:
    #  sims/exampleScript  -> <empty>
    #  sims/exampleScript/ -> <empty>
    #  sims/example.exp/script_001  -> example.exp
    #  sims/a.fam/b.fam/script_001  -> a.fam/b.fam
    #  sims/a.fam/b.fam/script_001/ -> a.fam/b.fam
    # Note, this may be used below to fill in args.outfile
    args.expt_name = '/'.join(os.path.dirname(args.ens1).split('/')[1:])
    # readable version of expt_name
    args.expt_name_readable = args.expt_name if args.expt_name else 'Root'
    # get the Ensemble-set name from the first argument above
    if args.ens1.startswith('sims/'):
        args.sim_root = os.path.dirname(args.ens1)
        args.script_root = os.path.dirname(args.ens1).replace('sims', 'Scripts', 1)
    else:
        # not recognized -- give up
        args.sim_root = ''
        args.script_root = ''
    if args.verbose:
        print(f"{args.progname}: Inferring script base directory: `{args.script_root}'")

    # 3/ index
    # Note: it's OK if indexfile is given but the file is not found
    if args.index_off:
        print(f"{args.progname}: Not seeking an index file")
        # ensure over-ride
        args.indexfile = ''
    else:
        if not args.indexfile:
            # note, expt_name has sims/ removed already
            args.indexfile = os.path.join('Scripts', args.expt_name, 's_index.json')
            print(f"{args.progname}: Inferring index filename: `{args.indexfile}'")
        else:
            print(f"{args.progname}: Supplied index filename: `{args.indexfile}'")

    # the code for "make exp-reduce" (currently) supplies the index filename, which flows to here
    # bona fide experiments should have "s_index.json", but families will not.
    # expedient solution: don't insist on it, so that families can be reduced
    if args.indexfile:
        if os.access(args.indexfile, os.R_OK):
            print("%s: Index file `%s' was found." % (args.progname, args.indexfile))
        else:
            # we typically will look for an indexfile, but don't warn unless it's a .exp
            if args.expt_name_readable.endswith('.exp') or args.verbose:
                print("%s: Warning: Associated index file `%s' is not readable." % (args.progname, args.indexfile))
            print(f"{args.progname}: Proceeding without index file (not an error).")
            args.indexfile = None

    # 4/ outfile
    # the below enforces Sandbox conventions on output filenames
    if not args.outfile:
        args.outfile = os.path.join('sims', args.expt_name, 'reduce-%s.%s')
        print(f"{args.progname}: Inferring output file pattern: `{args.outfile}'")
    if args.outfile.count('%s') != 2:
        print(f"{args.progname}: Need two instances of '%s' in output file template")
        sys.exit(1)
    # ensure enclosing dir exists
    directory = os.path.dirname(args.outfile % ('dummy', 'txt'))
    if not os.path.exists(directory):
        os.makedirs(directory)

    # Finally!
    print(f"{args.progname}: Beginning computations for `{args.expt_name_readable}'")
    main(args)
    # measure elapsed time
    time_delta = time.time() - time_start
    print(f'{args.progname}: Done. Elapsed time: {time_delta:.3}s.')

    sys.exit(0)
