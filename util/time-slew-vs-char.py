#!/usr/bin/env python
r'''
time-slew-vs-char.py: reduce a pile of DRMs to slew-time vs. char-time

usage:
  time-slew-vs-char.py [ -O outfile ] [ -j N ] DRM [ ...]

where:
  DRM ... is a list of DRM pickles,
and:
  -O outfile gives a template (containing exactly two occurrences of %s) for
     file outputs.  This is optional - output will go to the sims/SCRIPT/...
     directory if not given.
  -j N means to use N parallel workers to process the files.  By default,
     about 2/3 of the available cores will be used (20 on aftac1, 30 on aftac2).
     If N = 0 or 1, no parallel workers are used, which is helpful for debugging.

Typical usage:
  util/time-slew-vs-char.py sims/HabEx_4m_TSDDtemp_top130DD_dmag26p0_20180408/drm/*

Note: This started as a one-off.  I adapted the DRM data reduction script because
I wanted to see if a functional "shell", suitable for one-offs, could be extracted 
from the rather complex full-reduction setup.

turmon jan 2018, oct 2018
'''

from __future__ import division
from __future__ import print_function
import argparse
import sys
import glob
import time
import json
import os
import gc
import csv
import warnings
from functools import partial
from collections import defaultdict, Counter
import cPickle as pickle
#import multiprocessing.dummy
import multiprocessing as mproc
import numpy as np
from scipy import stats
from scipy.interpolate import interp1d
import astropy.units as u

##
## Globals
##

# global verbosity mode, also usable for debugging print's
VERBOSITY = 0

# number of CPUs
N_CPU = mproc.cpu_count()

# filter incoming DRMs to characterization-only?  (disabled)
#MODE = 'char'
#MODE = 'det'
MODE = '*all*'

########################################
###
###  Utility Functions
###
########################################

def strip_units(x):
    r'''Strip astropy units from x.'''
    # TODO: allow coercing units to a supplied value
    if hasattr(x, 'value'):
        return x.value
    else:
        return x

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

class SimulationRun(object):
    r'''Load and summarize a simulation: one DRM and its corresponding SPC.'''
    def __init__(self, f):
        # allow creating a dummy object so that its properties may be queried
        if f is None:
            self.drm = []
            self.spc = defaultdict(list)
            self.name = 'dummy'
            self.Nstar = 0
            return
        # disabling gc during object construction speeds up by ~30% (12/2017, py 2.7.14)
        gc.disable()
        drm = pickle.loads(open(f).read())
        gc.enable()
        # sometimes, skip some drms - generally unused.
        #if args.drm1 and len(drm) > 1: continue
        #if args.drm2 and len(drm) < 2: continue
        # allow to filter the events - disabled at present
        if 'char' in MODE:
            drm_filter = [d for d in drm if 'char_mode' in d]
            drm = drm_filter
        elif 'det' in MODE:
            drm_filter = [d for d in drm if 'det_status' in d]
            drm = drm_filter
        # load a spc file
        g = f.replace('pkl', 'spc').replace('/drm/', '/spc/')
        if os.path.isfile(g):
            spc = pickle.loads(open(g).read())
        else:
            raise ValueError('Could not find a .spc file to match DRM <%s>' % f)
        # set up object state
        self.name = f
        self.seed = int(os.path.splitext(os.path.basename(f))[0])
        self.Nstar = len(spc['Name'])
        self.spc = spc
        self.drm = drm
        self.summary = None # place-holder
    
    def extract_char_slew(self):
        r'''Extract char and slew times from DRM.

        Algorithm:
        For each event class, keep a count of each event of that
        type across the DRM.  That count will *later* be made into 
        a histogram.  In contrast to other routines, we return a set of scalars.
        '''

        # keep track of all event counts as a list
        char_time = []
        slew_time = []
        obs_num = []
        # DRM-FMT
        for obs in self.drm:
            # Process a detection
            if 'det_time' in obs or 'det_info' in obs:
                pass
            # Process a characterization
            if 'char_mode' in obs or 'char_info' in obs:
                if 'char_info' in obs:
                    char_info_1 = obs['char_info'][0]
                else:
                    char_info_1 = obs
                char_time_1 = strip_units(char_info_1['char_time']) # in days
                slew_time_1 = strip_units(char_info_1.get('slew_time', 0.0)) # may not exist
                # skip "time = 0" chars: they are an artifact
                # note: if slew_time was not there, it will have been set to 0 above
                if char_time_1 > 0:
                    char_time.append(char_time_1)
                    slew_time.append(slew_time_1)
                    obs_num.append(len(slew_time))
        # return value
        rv = dict(
            char_time=char_time,
            slew_time=slew_time,
            obs_num=obs_num,
            seed=[self.seed]*len(char_time))
        return rv

    def summarize(self, econo=True):
        r'''Find the summary of the sim as a dictionary held within the object.

        The convention is that the summary is built up by calling a series of analytic
        routines, each of which returns a dictionary of summary information.  The overall
        summary dictionary is a union of each individual summary.
        If econo, delete the DRM and keep only the summary.'''
        # this dict holds reductions for the current sim
        summary = {}
        # fold in event counts
        summary.update(self.extract_char_slew())
        # delete the base data if asked
        if econo:
            self.drm = None
            self.spc = None
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
    sim = SimulationRun(f)
    return sim.summarize()

class EnsembleSummary(object):
    r'''Compute, store, and dump summary information for an ensemble of many simulations.'''
    def __init__(self, in_files, args, lazy=True):
        '''Load an ensemble of simulations.'''
        # filter the .spc files out: some dirs contain both
        sim_files = [f for f in in_files if not f.endswith('.spc')]
        # save some useful state
        self.args = args
        self.sim_files = sim_files
        self.Ndrm_actual = len(sim_files) # = 0 if no sims
        # load on __init__ vs. load on reduce
        if lazy:
            # lazy => load the sims only when doing the reduce()
            # just get one run, for the SPC
            if len(sim_files) > 0:
                sims = [SimulationRun(sim_files[0])]
        else:
            # Load the sims all at once
            # THIS IS NO LONGER THE RECOMMENDED PATH
            sims = map(SimulationRun, sim_files)
        # Save the sims in the object
        # As a convenience, we insert a "phony" empty sim (DRM/SPC) even if
        # there were no input DRMs.  This allows all the histograms to be
        # initialized to their correct length, etc.
        # NB: Ndrm_actual can = 0, but self.sims will still contain 1 element
        if len(sim_files) == 0:
            self.sims = [SimulationRun(None)]
        else:
            self.sims = sims

    def load_and_reduce(self):
        r'''Load sim drm/spc, reduce each sim, accumulate summaries across sims.

        Each dict gives one summary statistic over that single sim.'''

        # In general, creates a pool of workers (separate unix processes)
        # (but: if jobs <= 1, it uses ordinary python map() and does no multiprocessing)
        with WorkerMap(self.args.jobs) as map_function:
            # map the load-and-reduce function over each file
            # reductions is a list of dicts containing summaries
            reductions = map_function(partial(outer_load_and_reduce, verb=self.args.verbose),
                                      self.sim_files)
        # hacky fix for no-drm case
        if len(self.sim_files) == 0:
            reductions = outer_load_and_reduce(None)
        # re-group the above reductions across sims
        # result is a dict containing reduced data, stored as self.summary
        self.regroup_and_accum(reductions)

    def regroup_and_accum(self, reductions):
        r'''Accumulate various summaries across the ensemble.

        Nothing is returned: result is placed in the object state.'''
        
        # list of attributes to accumulate
        attrs = [
            'seed', 'obs_num', 'char_time', 'slew_time', 
            ]
        # flatten the reductions from [drm][attribute] to [attribute]
        # summary is a dictionary of lists
        summary = {}
        for attr in attrs:
            summary[attr] = []
            for r in reductions:
                summary[attr].extend(r[attr])
        # record these summaries in the object
        self.summary = summary

    def dump(self, args):
        r'''Dump reduced data to files.'''

        def ensure_permissions(fn):
            r'''Ensure correct permissions on the named data file.  
            We use rw-rw-r-- = 664 (octal), to allow group-write.'''
            try:
                os.chmod(fn, 0o664)
            except OSError:
                pass # e.g., don't own the file

        fn = args.outfile % ('char-slew-time', 'csv')
        print('\tDumping to %s' % fn)
        # qoi = quantities-of-interest
        saved_fields = ['seed', 'obs_num', 'char_time',  'slew_time']
        Nrow = len(self.summary[saved_fields[0]])
        with open(fn, 'w') as csvfile:
            w = csv.DictWriter(csvfile, fieldnames=saved_fields)
            w.writeheader()
            for i in range(Nrow):
                # dictionary mapping field -> value -- everything is a scalar here
                d = {f:self.summary[f][i] for f in saved_fields}
                w.writerow(d)
        ensure_permissions(fn)
        

def main(args):
    print('%s: Loading %d file patterns.' % (args.progname, len(args.infile)))
    lazy = True
    ensemble = EnsembleSummary(args.infile, args, lazy=lazy)
    if ensemble.Ndrm_actual == 0:
        print('%s: Warning: no actual DRMs present.' % (args.progname, ))
    # save a reference to the universe (stars and their attributes)
    ensemble.spc = ensemble.sims[0].spc
    ensemble.Nstar = ensemble.sims[0].Nstar
    print('%s: Reducing.' % args.progname)
    ensemble.load_and_reduce()
    print('%s: Dumping.' % args.progname)
    ensemble.dump(args)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Summarize DRMs.",
                                     epilog='')
    parser.add_argument('drm', metavar='DRM', nargs='*', default='.',
                            help='drm file list')
    parser.add_argument('-O', '--outfile', type=str, default='',
                            help='Output file template.')
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
    
    # do it like this for now
    args.infile = args.drm
    # args.infile = 'drms/%s.pkl' % '*'

    # get the experiment name from the directory - this is brittle,
    # but have to do something
    args.expt_name = os.path.dirname(args.infile[0]).split('/')[-2]

    # best practice is to explicitly give outfile
    if not args.outfile:
        args.outfile = ('sims/%s/adhoc' % args.expt_name) + '-%s.%s'
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
