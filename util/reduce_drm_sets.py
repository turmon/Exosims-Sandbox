#!/usr/bin/env python
r'''
reduce_drm_sets.py: reduce a set of simulation-ensembles to summary CSV files

usage:
  reduce_drm_sets.py [ -O outfile ] [ -j N ] [ -i indexfile] ENS [...]

where:
  ENS ... is a list of simulation *directories* ("ensembles"),
and:
  -O outfile gives a template (containing exactly two occurrences of %s) for
     file outputs.  This is optional - output will go to the sims/SCRIPT/...
     directory if not given.
  -j N means to use N parallel workers to process the files.  By default,
     about 2/3 of the available cores will be used.
     If N = 0 or 1, no parallel workers are used, which is helpful for debugging.
  -i indexfile names a JSON index file that associates experiment names with
     parameter values for that experiment.  Giving this allows output of 
     summary information that is also labeled with the corresponding parameter values.

This program rolls up the summary information in ENS/reduce_info.csv, for each 
named ENS, and places the cumulative summary in a CSV file following the
given filename template. It also tabulates some yield metrics and places in related CSV files.

Typically a given ENS will be either a simulation ensemble (and ENS/drm will exist),
or ENS will be a family or experiment (without ENS/drm). 
In either case, the required information will be in ENS/reduce-info.csv.

Typical usage:
  util/reduce_drm_sets.py sims/HabEx_4m_TSDDtemp_top130DD_dmag26p0_20180408.exp/*

turmon nov 2018, feb 2022
'''

import argparse
import sys
import glob
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
    def __init__(self, f):
        # allow creating a dummy object so that its properties may be queried
        self.summary = None # place-holder for later
        if f is None:
            self.name = 'dummy'
            self.Nens = 0
            self.info = {}
            return
        info = self.read_sim_summary(f)
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
            ('experiment', str.lstrip, ''), 
            ('detections_earth_all', float, 0.0),
            ('detections_earth_unique', float, 0.0),
            ('chars_earth_unique', float, 0.0),
            ('chars_earth_strict', float, 0.0),
            # not used in the roll-up
            ('detections_unique_mean', float, 0.0),
            ('chars_unique_mean', float, 0.0), 
            ]
        props = {}
        for key, converter, nullval in prop_map:
            props[key] = converter(raw_props.get(key, nullval))
        return props

    def read_sim_summary(self, d):
        r'''Summarize one simulation directory (ensemble) into a dict.'''
        # grab the summary data for d
        info_fn = os.path.join(d, 'reduce-info.csv')
        try:
            with open(info_fn) as f:
                info_items = csv.DictReader(f);
                info = self.convert_sim_summary(next(info_items)) # it is a 1-line csv
                # if clause catches (and excludes) placeholder reduce-info's
                # that were created by the Dakota stub that generates scripts
                if len(info) < 2:
                    # print(f'\tSkipping: {info_fn}')
                    # this will later exclude the record
                    info = dict()
        except IOError:
            info = {}
        # Number of paths generated
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
def outer_load_and_reduce(f, verb=0):
    r'''Load a sim and summarize it into a dict.

    This must be present at the outer scope of the file so it can be loaded
    by a separate process that is created by the multiprocessing module.'''
    if verb > 0:
        print('Processing <%s> in pid #%d' % (f, os.getpid()))
    ensemble = EnsembleRun(f)
    return ensemble.summarize()

class EnsembleSummary(object):
    r'''Compute, store, and dump summary information for a set of ensembles.'''
    def __init__(self, in_files, args, lazy=True):
        '''Load an ensemble of simulations.'''
        # load index file, if given
        self.load_index(args)
        # filter: only retain input directories with reduced DRM-sets
        ens_files = [f for f in in_files if os.path.isfile(os.path.join(f, 'reduce-info.csv'))]
        # save some useful state
        self.args = args
        self.ens_files = ens_files
        self.Nens = len(ens_files) # = 0 if no ensembles

    def load_index(self, args):
        r'''Load index mapping experiment name to parameter values.'''
        # load index file, if given
        if args.indexfile:
            print('%s: Loading index file.' % args.progname)
            with open(args.indexfile, 'r') as fp:
                index = json.load(fp)
        else:
            index = []
        # copy index into a dict-of-dicts, EXP_NAME -> {param1:value1, param2:value2, ...}
        self.index = dict()
        for s in index:
            s1 = {k:v for k,v in s.items() if not k.endswith('_name')}
            self.index[s['run_name']] = s1
        # also create a scalarized version of the index, for CSV output
        # this is the same mapping as "index", but it expands value-lists into sequences of named scalars
        index_csv = dict()
        for exp_name, d in self.index.items():
            index_csv[exp_name] = OrderedDict()
            for param, value in d.items():
                if isinstance(value, list):
                    # expand the vector parameter "value" element-by-element
                    for inx, v1 in enumerate(value):
                        index_csv[exp_name]['%s%d' % (param, inx+1)] = v1
                else:
                    index_csv[exp_name][param] = value
        # save it in the object
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
            reductions = list(map_function(partial(outer_load_and_reduce, verb=self.args.verbose),
                                      self.ens_files))
        # hacky fix for no-drm case
        # note: n_valid <= len(self.ens_files)
        n_valid = len([1 for r in reductions if len(r) > 0])
        if n_valid == 0:
            reductions = outer_load_and_reduce(None)
        # need to save the original data for full-results tabular output
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

        # 2: take means, std's, quantiles of the various attributes
        #    some QOIs can have NaNs (eg, char_snr => from empty RpL bins;
        #    h_star_det_earth_frac => from stars without Earths)
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
        # some additional info
        extra_info = dict(
            user=os.environ['USER'],
            runtime=time.strftime("%Y-%m-%d_%H:%M"),
            simtime=simtimes[-1] if simtimes else '2000-01-01_00:00',
            experiment=args.expt_name_readable,
            experiment_size=len(reductions), # only already-reduced ensembles
            )
        summary.update(extra_info)
        # record this cross-ensemble summary in the object
        self.summary = summary


    def dump_results_worker(self, args, extension, first_fields):
        r'''Dump reduced data to files.'''

        fn = args.outfile % (extension, 'csv')
        print('\tDumping to %s' % fn)
        # list of field-name in order they should be dumped
        saved_fields = [
            'chars_earth_unique',
            'detections_earth_all',
            'detections_earth_unique',
            'ensemble_size',
            ]
        # take parameter-list, if any, from first param entry in self.index_csv
        has_index = len(self.index_csv) > 0
        param_fields = list(self.index_csv[list(self.index_csv.keys())[0]].keys()) if has_index else []
        # fields to dump includes both
        all_fields = first_fields + param_fields + saved_fields
        # TODO: Extend field map with varying parameters
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
                    exp_name = r['experiment']
                    try:
                        exp_params = self.index_csv[exp_name]
                    except KeyError:
                        # happens if there are ensembles in the .exp dir, that are not
                        # in s_index.json
                        print(f'{args.progname}: params for {exp_name} not in s_index.json, is it up-to-date?')
                        exp_params = defaultdict(lambda: "?")
                    for f in param_fields:
                        d[f] = exp_params[f]
                w.writerow(d)
        ensure_permissions(fn)
        
    def dump_results(self, args):
        self.dump_results_worker(args, 'yield', [])
        self.dump_results_worker(args, 'yield-plus', ['experiment'])

    def dump_summary(self, args):
        r'''Dump overall, one-line summary data to a file.'''

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


def main(args):
    print(f'{args.progname}: Examining {len(args.infile)} file patterns.')
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
    print(f'{args.progname}: Done.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Summarize Experiments = Ensemble Sets.",
                                     epilog='')
    parser.add_argument('ens', metavar='ENS', nargs='*', default=[],
                            help='ensemble list')
    parser.add_argument('-O', '--outfile', type=str, default='',
                            help='Output file template.')
    parser.add_argument('-i', '--indexfile', type=str, default='',
                            help='Index file name.')
    parser.add_argument('-v', help='verbosity', action='count', dest='verbose',
                            default=0)
    parser.add_argument('-j', '--jobs', help=' # parallel jobs, default = %(default)d',
                      type=int, dest='jobs', default=int(N_CPU*0.65))
    args = parser.parse_args()
    
    # set umask in hopes that files will be group-writable
    os.umask(0o002)

    VERBOSITY = args.verbose

    # program name, for convenience
    args.progname = os.path.basename(sys.argv[0])
    
    args.infile = args.ens
    if len(args.infile) == 0:
        print('%s: Error: Need at least one input ensemble.' % args.progname)
        sys.exit(1)

    # the code for "make exp-reduce" supplies the index filename, which flows to here
    # bona fide experiments should have the s_index.json file, but families will not.
    # expedient solution: don't insist on it, so that families can be reduced
    if args.indexfile and not os.access(args.indexfile, os.R_OK):
        print("%s: Warning: Associated index file `%s' is not readable." % (args.progname, args.indexfile))
        print("%s: Proceeding without index file." % (args.progname, ))
        args.indexfile = None

    # get the experiment name from the directory - brittle, but expedient.
    # examples of operation:
    #  sims/exampleScript  -> <empty>
    #  sims/exampleScript/ -> <empty>
    #  sims/example.exp/script_001  -> example.exp
    #  sims/a.fam/b.fam/script_001  -> a.fam/b.fam
    #  sims/a.fam/b.fam/script_001/ -> a.fam/b.fam
    # Note, this may be used below to fill in args.outfile
    infile0 = args.infile[0]
    if infile0.endswith('/'):
        infile0 = infile0[:-1] # kill trailing / if supplied
    args.expt_name = '/'.join(os.path.dirname(infile0).split('/')[1:])

    # readable version of expt_name
    args.expt_name_readable = args.expt_name if args.expt_name else 'Root'

    # get the Ensemble-set name from the first argument above
    args.script_root = os.path.dirname(infile0).replace('sims/', 'Scripts/')
    # print('%s: Inferring script file: %s' % (args.progname, args.script_root))

    # best practice is to explicitly give outfile, this is the backup
    if not args.outfile:
        args.outfile = os.path.join('sims', args.expt_name, 'reduce-%s.%s')
        print("%s: Inferring output file pattern: %s." % (args.progname, args.outfile))
    if args.outfile.count('%s') != 2:
        print('%s: Need two instances of %%s in output file template' % args.progname)
        sys.exit(1)
    # ensure enclosing dir exists
    directory = os.path.dirname(args.outfile % ('dummy', 'txt'))
    if not os.path.exists(directory):
        os.makedirs(directory)

    infile_print = args.infile[0] + (' ...' if len(args.infile) > 1 else '')
    print('%s: Reducing %s to %s' % (args.progname, infile_print, args.outfile))

    main(args)
    sys.exit(0)
