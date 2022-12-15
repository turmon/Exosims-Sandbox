#!/usr/bin/env python
r'''
drm-grep.py: search through a pile of DRMs for named fields

usage:
  drm-grep.py [-1snv] [ -m MATCH ] [-a ATTR] [ -j N ] DRM ...

where:
  DRM ... is a list of DRM pickles,
and:
  -a ATTR means: output the keyword ATTR (repeat -a for multiple keywords,
     or separate their names by commas)
  -S ATTR means: output the SPC attribute ATTR for obs[star_ind]
  -P ATTR means: output the SPC attribute ATTR for obs[plan_inds[0]]
  -m MATCH means: only produce output if MATCH is a keyword in the observation
  -1 means: output the CSV header ("line 1")
  -s means: output the DRM seed number
  -n means: output the DRM observation number

or the less-useful options:
  -v increases verbosity (output is to stderr)
  -j N means to use N parallel workers to process the files.  By default,
     about 2/3 of the available cores will be used (20 on aftac1, 30 on aftac2).
     If N = 0 or 1, no parallel workers are used, which is helpful for debugging.

Extension:
Some attributes are nested.  To get at nested attributes, such as:
  drm[0]['char_mode']['lam']
use:
  -a char_mode.lam
This extends to multiply-nested attributes, and to lists.  The initial plan_ind is
  -a plan_inds.0
The last phase angle is:
  -a char_params.phi.-1

Typical usage:
  # (all arrival_time's)
  util/drm-grep.py -a arrival_time sims/HabEx_4m_TSDDtemp_top130DD_dmag26p0_20180408/drm/*

  # (char_time, but only looking at characterizations)
  util/drm-grep.py -m char_time -a char_time sims/HabEx_4m_TSDDtemp_top130DD_dmag26p0_20180408/drm/*

turmon apr 2019
'''

#####
#####
#####  THIS FILE HAS BEEN SUPERSEDED BY drm-tabulate.py
#####
#####

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
import six.moves.cPickle as pickle
#import multiprocessing.dummy
import multiprocessing as mproc
import numpy as np
import astropy.units as u
from six.moves import map
from six.moves import range

##
## Globals
##

# global verbosity mode, also usable for debugging print's
VERBOSITY = 0

# number of CPUs
N_CPU = mproc.cpu_count()

# unpickling python2/numpy pickles within python3 requires this
PICKLE_ARGS = {} if sys.version_info.major < 3 else {'encoding': 'latin1'}

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
            # map function is the python map() builtin (forced to materialize the list)
            self.map_function = lambda f,x: list(map(f,x))
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
    def __init__(self, f, args):
        # allow creating a dummy object so that its properties may be queried
        if f is None:
            self.drm = []
            self.spc = defaultdict(list)
            self.name = 'dummy'
            self.Nstar = 0
            return
        # disabling gc during object construction speeds up by ~30% (12/2017, py 2.7.14)
        gc.disable()
        drm = pickle.load(open(f, 'rb'), **PICKLE_ARGS)
        gc.enable()
        # load a spc file - only if needed to lookup -S, -P attrs
        if args.load_spc:
            g = f.replace('pkl', 'spc').replace('/drm/', '/spc/')
            if os.path.isfile(g):
                spc = pickle.load(open(g, 'rb'), **PICKLE_ARGS)
            else:
                raise ValueError('Could not find a .spc file to match DRM <%s>' % f)
        else:
            spc = None
        # set up object state
        self.name = f
        self.seed = int(os.path.splitext(os.path.basename(f))[0])
        # self.Nstar = len(spc['Name'])
        self.spc = spc # is None, if SPC not needed
        self.drm = drm
        self.summary = None # place-holder
    
    def extract_value(self, obs, attr):
        r'''Extract a value, attr, from an observation, obs, allowing recursive lookup.'''
        if '.' not in attr:
            # if we can convert the attribute to an int, do so, and it will
            # be used as a list index, not a dict lookup
            try:
                attr_use = int(attr)
            except:
                attr_use = attr
            # look up attr in obs
            try:
                rv = strip_units(obs[attr_use])
            except:
                rv = None
            return rv
        else:
            segments = attr.split('.')
            # if we can convert the first segment to an int, do so, and it will
            # be used as a list index, not a dict lookup
            try:
                seg0 = int(segments[0])
            except:
                seg0 = segments[0]
            # recursive extract
            return self.extract_value(obs[seg0], '.'.join(segments[1:]))

    def extract_attrs(self, match, attrs, sattrs, pattrs):
        r'''Extract attributes from DRM.

        Algorithm:
        For each event class, keep a count of each event of that
        type across the DRM.
        '''

        # keep track of all event counts as a list
        attr_vals = {attr: [] for attr in (attrs + sattrs + pattrs)}
        obs_num = []
        # DRM-FMT
        for obs in self.drm:
            if match and match not in obs:
                # skip this observation
                continue
            for attr in attrs:
                attr_val_1 = self.extract_value(obs, attr)
                attr_vals[attr].append(attr_val_1)
            for attr in sattrs:
                sind = obs['star_ind']
                attr_val_1 = strip_units(self.spc[attr][sind])
                attr_vals[attr].append(attr_val_1)
            for attr in pattrs:
                plan_ind = obs['plan_inds'][0]
                attr_val_1 = strip_units(self.spc[attr][plan_ind])
                attr_vals[attr].append(attr_val_1)
            # obtain this for possible use later
            if 'ObsNum' in obs:
                obs_num_1 = obs['ObsNum']
            else:
                obs_num_1 = obs['Obs#']
            obs_num.append(obs_num_1)
        # return value -- fixed fields
        rv = dict(
            obs_num=obs_num,
            seed=[self.seed]*len(attr_vals[attrs[0]]))
        # add the dynamic fields
        # NB: we could instead tensorize over planets here (?)
        for attr in (attrs + sattrs + pattrs):
            rv[attr] = attr_vals[attr]
        return rv


    def summarize(self, args, econo=True):
        r'''Find the summary of the sim as a dictionary held within the object.

        The convention is that the summary is built up by calling a series of analytic
        routines, each of which returns a dictionary of summary information.  The overall
        summary dictionary is a union of each individual summary.
        If econo, delete the DRM and keep only the summary.'''
        # this dict holds attribute-lists for the current sim
        summary = self.extract_attrs(args.match, args.attrs, args.sattrs, args.pattrs)
        # delete the base data if asked
        if econo:
            self.drm = None
        # keep a reference in the object
        self.summary = summary
        # also return the summary-dictionary
        return summary


# NB: Present in outer scope as a helper for load_and_reduce() below.
def outer_load_and_reduce(f, verb=0, args=None):
    r'''Load a sim and summarize it into a dict.

    This must be present at the outer scope of the file so it can be loaded
    by a separate process that is created by the multiprocessing module.'''
    if verb > 0:
        print('Processing <%s> in pid #%d' % (f, os.getpid()))
    sim = SimulationRun(f, args)
    return sim.summarize(args)

class EnsembleSummary(object):
    r'''Compute, store, and dump summary information for an ensemble of many simulations.'''
    def __init__(self, in_files, args):
        '''Prepare to load an ensemble of simulations.'''
        # filter the .spc files out: some dirs contain both
        sim_files = [f for f in in_files if not f.endswith('.spc')]
        # save some useful state
        self.args = args
        self.sim_files = sim_files
        self.Ndrm_actual = len(sim_files) # = 0 if no sims

    def load_and_reduce(self):
        r'''Load sim drm/spc, reduce each sim, accumulate summaries across sims.

        Each dict gives one summary statistic over that single sim.'''

        # In general, creates a pool of workers (separate unix processes)
        # (but: if jobs <= 1, it uses ordinary python map() and does no multiprocessing)
        with WorkerMap(self.args.jobs) as map_function:
            # map the load-and-reduce function over each file
            # reductions is a list of dicts containing summaries
            reductions = map_function(partial(outer_load_and_reduce, args=args, verb=self.args.verbose),
                                      self.sim_files)
        
        # hacky fix for no-drm case
        if len(self.sim_files) == 0:
            reductions = outer_load_and_reduce(None)
        # re-group the above reductions across sims
        # result is a dict containing reduced data, stored as self.summary
        self.regroup_and_accum(reductions, (args.attrs + args.sattrs + args.pattrs))

    def regroup_and_accum(self, reductions, attrs):
        r'''Accumulate various summaries across the ensemble.

        Nothing is returned: result is placed in the object state.'''
        
        # list of attributes to accumulate
        attrs_used = [
            'seed', 'obs_num', 
            ]
        attrs_used.extend(attrs)
        # flatten the reductions from [drm][attribute] to [attribute]
        # summary is a dictionary of lists
        summary = {}
        for attr in attrs_used:
            summary[attr] = []
            for r in reductions:
                summary[attr].extend(r[attr])
        # record these summaries in the object
        self.summary = summary

    def dump(self, args, outfile):
        r'''Dump reduced data to output.'''

        # qoi = quantities-of-interest
        saved_fields = []
        if args.seed:
            saved_fields.append('seed')
        if args.obs_num:
            saved_fields.append('obs_num')
        saved_fields.extend(args.attrs)
        saved_fields.extend(args.sattrs)
        saved_fields.extend(args.pattrs)
        Nrow = len(self.summary[saved_fields[0]])

        csv_opts = {}
        if args.delimiter:
            csv_opts['delimiter'] = args.delimiter
        w = csv.DictWriter(outfile, fieldnames=saved_fields, **csv_opts)
        if args.header:
            w.writeheader()
        for i in range(Nrow):
            # dictionary mapping field -> value -- everything is a scalar here
            d = {f:self.summary[f][i] for f in saved_fields}
            w.writerow(d)


def main(args, outfile):
    if VERBOSITY:
        sys.stderr.write('%s: Loading %d file patterns.\n' % (args.progname, len(args.infile)))
    ensemble = EnsembleSummary(args.infile, args)
    if ensemble.Ndrm_actual == 0:
        print('%s: Warning: no actual DRMs present.' % (args.progname, ))
    if VERBOSITY:
        sys.stderr.write('%s: Reducing.\n' % args.progname)
    ensemble.load_and_reduce()
    # print '%s: Dumping.' % args.progname
    # default is ~70, which means small-sized vectors are line-wrapped in the CSV
    np.set_printoptions(linewidth=1000)
    ensemble.dump(args, outfile)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Look for attributes in DRMs.",
                                     epilog='')
    parser.add_argument('drm', metavar='DRM', nargs='*', default='.',
                            help='drm file list')
    parser.add_argument('-m', '--match', type=str, default='',
                            help='Attribute to match.')
    parser.add_argument('-a', '--attrs', type=str, action='append', required=True, 
                            dest='attrs_orig', help='Attribute to extract, one for each -a.')
    parser.add_argument('-S', '--star', type=str, action='append', default=[],
                            dest='sattrs_orig', help='Star attribute to extract, one for each -S.')
    parser.add_argument('-P', '--plans', type=str, action='append', default=[],
                            dest='pattrs_orig', help='Planet attribute to extract, one for each -P.')
    parser.add_argument('-v', help='verbosity', action='count', dest='verbose',
                            default=0)
    parser.add_argument('-1', help='Supply header line', action='store_true', dest='header',
                            default=False)
    parser.add_argument('-d', '--delimiter', help='CSV field delimiter', type=str)
    parser.add_argument('-s', '--seed', help='Output the DRM seed', action='store_true')
    parser.add_argument('-n', '--obs_num', help='Output the observation number', action='store_true')
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

    # allow comma-separated attr/sattr/pattr lists
    attrs = []
    for attr in args.attrs_orig:
        attrs.extend(attr.split(','))
    args.attrs = attrs

    attrs = []
    for attr in args.sattrs_orig:
        attrs.extend(attr.split(','))
    args.sattrs = attrs

    attrs = []
    for attr in args.pattrs_orig:
        attrs.extend(attr.split(','))
    args.pattrs = attrs

    # do we need to load the SPC file?
    args.load_spc = args.sattrs or args.pattrs

    # get the experiment name from the directory - this is brittle,
    # but have to do something
    try:
        args.expt_name = os.path.dirname(args.infile[0]).split('/')[-2]
    except:
        args.expt_name = 'Simulation'

    outfile = sys.stdout

    main(args, outfile)
    sys.exit(0)
