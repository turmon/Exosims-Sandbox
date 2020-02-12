#!/usr/bin/env python
r'''
reduce_drm_sets.py: reduce a set of simulation-ensembles to a summary CSV file

usage:
  reduce_drm_sets.py [ -O outfile ] [ -j N ] [ -i indexfile] ENS [...]

where:
  ENS ... is a list of simulation *directories* ("ensembles"),
and:
  -O outfile gives a template (containing exactly two occurrences of %s) for
     file outputs.  This is optional - output will go to the sims/SCRIPT/...
     directory if not given.
  -j N means to use N parallel workers to process the files.  By default,
     about 2/3 of the available cores will be used (20 on aftac1, 30 on aftac2).
     If N = 0 or 1, no parallel workers are used, which is helpful for debugging.
  -i indexfile names a JSON index file that associates experiment names with
     parameter values for that experiment.  Giving this allows output of 
     summary information that is also labeled with the corresponding parameter values.

Note, ENS/drm should exist for each given ENS, and ENS/reduce_info.csv is 
read for each ENS.  If these directories/files do not exist, for example if
"make ... reduce" has not been run, then the ENS will be skipped.

Typical usage:
  util/reduce_drm_sets.py sims/HabEx_4m_TSDDtemp_top130DD_dmag26p0_20180408.exp/*

turmon nov 2018
'''

from __future__ import division
import argparse
import sys
import glob
import time
import os
import copy
import csv
import json
import warnings
from functools import partial
from collections import OrderedDict
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
        # properties we want to convert: key, converter
        # (we don't use all of these now)
        prop_map = [
            ('ensemble_size', int),
            ('runtime', str), # actually a date
            ('experiment', str), 
            ('detections_earth_all', float),
            ('detections_earth_unique', float),
            ('detections_unique_mean', float),
            ('chars_unique_mean', float),
            ('chars_earth_unique', float),
            ]
        props = {}
        for key, converter in prop_map:
            props[key] = converter(raw_props[key])
        return props

    def read_sim_summary(self, d):
        r'''Summarize one simulation directory (ensemble) into a dict.'''
        path_count = len(glob.glob(os.path.join(d, 'path/[0-9]*-cume/')))
        # grab the summary data for d
        info_fn = os.path.join(d, 'reduce-info.csv')
        try:
            with open(info_fn) as f:
                info_items = csv.DictReader(f);
                info = self.convert_sim_summary(info_items.next()) # it is a 1-line csv
            info['path_count'] = path_count
        except IOError:
            info = {}
        # Number of paths generated
        return info

    def extract_info(self):
        r'''Extract summary info about the Ensemble.'''
        return copy.deepcopy(self.info)

    def summarize(self, econo=True):
        r'''Find the summary of the Ensemble as a dictionary held within the object.

        The convention is that the summary is built up by calling a series of analytic
        routines, each of which returns a dictionary of summary information.  The overall
        summary dictionary is a union of each individual summary.
        If econo, delete the info and keep only the summary.'''
        # this dict holds reductions for the current sim
        summary = {}
        # fold in info -- this is trivial at the moment, pending more extensive summary
        summary.update(self.extract_info())
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
        print 'Processing <%s> in pid #%d' % (f, os.getpid())
    ensemble = EnsembleRun(f)
    return ensemble.summarize()

class EnsembleSummary(object):
    r'''Compute, store, and dump summary information for a set of ensembles.'''
    def __init__(self, in_files, args, lazy=True):
        '''Load an ensemble of simulations.'''
        # load index file, if given
        self.load_index(args)
        # filter: only retain directories with DRM-sets
        ens_files = [f for f in in_files if os.path.isdir(os.path.join(f, 'drm'))]
        # save some useful state
        self.args = args
        self.ens_files = ens_files
        self.Nens = len(ens_files) # = 0 if no ensembles

    def load_index(self, args):
        r'''Load index mapping experiment name to parameter values.'''
        # load index file, if given
        if args.indexfile:
            print '%s: Loading index file.' % args.progname
            with open(args.indexfile, 'r') as fp:
                index = json.load(fp)
        else:
            index = []
        # copy index into a dict-of-dicts, EXP_NAME -> {param1:value1, param2:value2, ...}
        self.index = dict()
        for s in index:
            s1 = {k:v for k,v in s.iteritems() if not k.endswith('_name')}
            self.index[s['run_name']] = s1
        # also create a scalarized version of the index, for CSV output
        # this is the same mapping as "index", but it expands value-lists into sequences of named scalars
        index_csv = dict()
        for exp_name, d in self.index.iteritems():
            index_csv[exp_name] = OrderedDict()
            for param, value in d.iteritems():
                if isinstance(value, list):
                    # expand the vector parameter "value" element-by-element
                    for inx, v1 in enumerate(value):
                        index_csv[exp_name]['%s%d' % (param, inx+1)] = v1
                else:
                    index_csv[k][param] = value
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
            reductions = map_function(partial(outer_load_and_reduce, verb=self.args.verbose),
                                      self.ens_files)
        # hacky fix for no-drm case
        if len(self.ens_files) == 0:
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
            'chars_unique_mean', 'chars_earth_unique',
            'ensemble_size', 'path_count',
            ]

        # 1: flatten the reductions from [ens][attribute] to [attribute]
        # accum is a dictionary of lists
        accum = {}
        for attr in attrs:
            accum[attr] = []
            for r in reductions:
                accum[attr].append(r[attr])

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
            experiment=args.expt_name,
            experiment_size=len(reductions), # only already-reduced ensembles
            )
        summary.update(extra_info)
        # record this cross-ensemble summary in the object
        self.summary = summary


    def dump_results_worker(self, args, extension, first_fields):
        r'''Dump reduced data to files.'''

        fn = args.outfile % (extension, 'csv')
        print '\tDumping to %s' % fn
        # list of field-name in order they should be dumped
        saved_fields = [
            'chars_earth_unique',
            'detections_earth_all',
            'detections_earth_unique',
            'ensemble_size',
            ]
        # take parameter-list, if any, from first param entry in self.index_csv
        has_index = len(self.index_csv) > 0
        param_fields = self.index_csv[self.index_csv.keys()[0]].keys() if has_index else []
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
                    exp_params = self.index_csv[exp_name]
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
        print '\tDumping to %s' % fn
        # list of (field-name-as-output, field-name-here),
        # in order they should be dumped
        saved_field_map = [
            ('user', 'user'),
            ('runtime', 'runtime'),
            ('experiment_size', 'experiment_size'),
            ('chars_earth_unique', 'chars_earth_unique_max'),
            ('detections_earth_all', 'detections_earth_all_max'),
            ('detections_earth_unique', 'detections_earth_unique_max'),
            ('path_count', 'path_count_sum'),
            ('ensemble_size', 'ensemble_size_sum'),
            ]
        with open(fn, 'w') as csvfile:
            w = csv.DictWriter(csvfile, fieldnames=[f for f,_ in saved_field_map])
            w.writeheader()
            # dictionary mapping field -> value -- everything is a scalar here
            d = {f_out:self.summary[f_have] for f_out, f_have in saved_field_map}
            w.writerow(d)
        ensure_permissions(fn)
        

def main(args):
    print '%s: Loading %d file patterns.' % (args.progname, len(args.infile))
    ensemble_set = EnsembleSummary(args.infile, args)
    if ensemble_set.Nens == 0:
        print '%s: Warning: no actual Ensembles present.' % (args.progname, )
    else:
        print '%s: Found %d ensembles.' % (args.progname, ensemble_set.Nens)
    # save a reference to the universe (stars and their attributes)
    print '%s: Reducing.' % args.progname
    ensemble_set.load_and_reduce()
    print '%s: Dumping result tabulation.' % args.progname
    ensemble_set.dump_results(args)
    print '%s: Dumping capsule summary.' % args.progname
    ensemble_set.dump_summary(args)


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
        print '%s: Error: Need at least one input ensemble.' % args.progname
        sys.exit(1)

    if args.indexfile and not os.access(args.indexfile, os.R_OK):
        print "%s: Error: Supplied index file `%s' is not readable." % (args.progname, args.indexfile)
        sys.exit(1)

    # get the experiment name from the directory - this is brittle,
    # but have to do something
    args.expt_name = os.path.dirname(args.infile[0]).split('/')[-1]

    # best practice is to explicitly give outfile
    if not args.outfile:
        args.outfile = ('sims/%s/reduce' % args.expt_name) + '-%s.%s'
    if args.outfile.count('%s') != 2:
        print '%s: Need two instances of %%s in output file template' % args.progname
        sys.exit(1)
    # ensure enclosing dir exists
    directory = os.path.dirname(args.outfile % ('dummy', 'txt'))
    if not os.path.exists(directory):
        os.makedirs(directory)

    infile_print = args.infile[0] + (' ...' if len(args.infile) > 1 else '')
    print '%s: Reducing %s to %s' % (args.progname, infile_print, args.outfile)

    main(args)
    sys.exit(0)
