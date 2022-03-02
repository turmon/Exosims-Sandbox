#!/usr/bin/env python
r'''
drm-tabulate.py: extract named fields from a pile of DRMs

usage:
  drm-tabulate.py [-1sn] [--json] [ -m ATTR ] [-a ATTR | -f FILE] DRM ...

where:
  DRM ... is a list of DRM pickles
and:
  -a ATTR: output the keyword ATTR (repeat -a OK, also see below)
  -f FILE: take ATTRs from FILE *instead of* -a options (see below)
  -m ATTR: only produce output if ATTR is a keyword in the observation
  -s: output the DRM seed number (__seed__, boolean)
  -n: output the DRM observation number (__obs_num__, boolean)
  -1: print line 1, the CSV header (__header__, boolean)

and these less-useful options:
  --json: output is JSON, rather than standard CSV
  -v: increase verbosity (output is to stderr)
  -j N: use N parallel workers (default, ~2/3 of cores)
        If N = 0 or 1, no parallelism: needed for debugging.

The main degree of freedom is attribute specification -- which is done
using either of two notations. Below, suppose "obs" is one entry in the DRM.

(1) Directly named attributes
    1a) Get a field in obs by giving: "-a arrival_time"
    1b) Drill into nested attributes with "."; for example
          obs['char_mode']['lam']
        is extracted with:
          -a char_mode.lam
        For lists, give the index number; obs['plan_inds'][0] is
          -a plan_inds.0 
        The last phase angle, obs['char_params'][-1], is:
          -a char_params.phi.-1
    1c) Supply a comma-separated list of such DRM fields in one -a
(2) Evaluated expressions
    2a) A python expression can be given, which is evaluated in the context 
        of variables named for the fields in "obs". To output a count of 
        detections using the obs['det_status'] list, use:
          -a "np.sum(det_status == 1)"

The resulting column can be custom-named with a label:attr construct,
          -a "lambda:char_mode.lam"
          -a "det_count:np.sum(det_status == 1)"
otherwise a basic generated name is used.
(But: Attributes in comma-separated expressions cannot be custom-named.)

These expressions can become complex, so the "-a ATTR" constructs can be 
placed in a JSON file and specified with -f FILE:

{
  "__match__": "char_info",
  "lambda": "char_mode.lam",
  "char_count": "np.sum(char_info[0]['char_status'] == 1)"
}

Recall JSON uses double-quotes, and sub-attributes are extracted with ...['attr'].
The key name (e.g., "char_count" above) is used as the column name for the CSV.
Above, "__match__" abbreviates the -m construct, so the JSON file can be
self-contained. Other program flags (like -s) can be given by Booleans in the 
JSON file as well. See the top of this usage note for the __attribute_name__
controlling each flag.

Typical usage:
  # each arrival_time
  util/drm-tabulate.py -a arrival_time sims/HabEx_4m_dmag26/drm/*

  # each char_time, looking only at chars
  util/drm-tabulate.py -m char_time -a char_time sims/HabEx_4m_dmag26/drm/*

  # number of successful chars, only for chars
  util/drm-tabulate.py -m char_info -a "np.sum(char_info[0]['char_status'] == 1)" sims/HabEx_4m_dmag26/drm/*

turmon apr 2019, feb 2022
'''


from __future__ import division
from __future__ import print_function
import argparse
import sys
import time
import json
import os
import re
import gc
import csv
import warnings
import collections
from functools import partial
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

class NumpyEncoder(json.JSONEncoder):
    r"""Custom JSON encoder for numpy types."""
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, u.quantity.Quantity):
            # note: it is possible to have a numpy ndarray wrapped in a Quantity,
            # and obj will be both a Quantity and an ndarray
            # for the moment, placing this ahead of ndarray works, although
            # Quantity vectors might complicate this.
            return obj.value
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, Time):
            # astropy Time -> time string
            return obj.fits # isot also makes sense here
        return json.JSONEncoder.default(self, obj)


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
            self.spc = collections.defaultdict(list)
            self.name = 'dummy'
            self.Nstar = 0
            return
        # disabling gc during object construction speeds up by ~30% (12/2017, py 2.7.14)
        gc.disable()
        drm = pickle.load(open(f, 'rb'), **PICKLE_ARGS)
        gc.enable()
        # load a spc file - disabled for now
        if False:
            g = f.replace('pkl', 'spc').replace('/drm/', '/spc/')
            if os.path.isfile(g):
                spc = pickle.load(open(g, 'rb'), **PICKLE_ARGS)
            else:
                raise ValueError('Could not find a .spc file to match DRM <%s>' % f)
        # set up object state
        self.name = f
        self.seed = int(os.path.splitext(os.path.basename(f))[0])
        # self.Nstar = len(spc['Name'])
        # self.spc = spc
        self.drm = drm
        self.summary = None # place-holder
    

    def extract_value(self, obs, attr, just_peeking=False):
        r'''Extract a value, attr, from an observation, obs, allowing recursive lookup.'''
        if not attr:
            # original obs[attr] is indeed in obs
            if just_peeking:
                return True
            # convert an astropy.quantity to a number, if possible
            rv = strip_units(obs)
            return rv
        else:
            # look up attr in obj
            segments = attr.split('.')
            # if segments[0] converts to int, use as list index, else, dict lookup
            try:
                seg0 = int(segments[0])
            except ValueError:
                seg0 = segments[0]
            # easiest to attempt the extraction here to see if it's OK
            try:
                _ = obs[seg0]
            except (KeyError, TypeError, IndexError):
                # KeyError: key not found
                # IndexError: list index bad
                # TypeError: list/dict confusion
                if just_peeking:
                    return False
                # announce the cause now, exception triggered below
                sys.stderr.write('Fatal: Value extraction failed while trying {}.\n'.format(
                    'list indexing' if type(seg0) is int else 'dict lookup'))
                sys.stderr.write('No "{}" found in obs object (or sub-object).\n'.format(seg0))
                sys.stderr.write('Object or sub-object is a {}:\n{}\n'.format(type(obs), repr(obs)))
            # extract from sub-object
            return self.extract_value(obs[seg0], '.'.join(segments[1:]), just_peeking=just_peeking)

    def extract_pseudo(self, obs, nobs, attr):
        r'''Extract pseudo-attribute, which can have special naming conventions.'''
        if attr == 'obs_num':
            if 'ObsNum' in obs:
                return obs['ObsNum']
            elif 'Obs#' in obs:
                return obs['Obs#']
            else:
                return nobs + 1
        elif attr == 'seed':
            return self.seed
        else:
            assert False, 'Should not be reached (attr = {})'.format(attr)
            
    def perform_eval(self, obs, attr):
        r'''Extract attribute by evaluating attr in the context of obs.'''
        try:
            val = eval(attr, {"np": np, **obs})
        except:
            sys.stderr.write('Error: Failed to eval() the given attribute.\n')
            sys.stderr.write('Evaluated string is surrounded by || below:\n|{}|\n'.format(attr))
            sys.stderr.write('Context is:\n')
            sys.stderr.write(repr(obs) + '\n')
            raise
        return val

    def extract_attrs(self, match, attrs, attr_names, attr_types):
        r'''Extract attributes from each obs in the DRM.

        Returns a dictionary mapping attributes to lists, one list entry
        per DRM observation.
        '''
        # accumulate values in these bins
        vals = {n: [] for n in attr_names}
        for nobs, obs in enumerate(self.drm):
            if match and not self.extract_value(obs, match, just_peeking=True):
                continue # no obs[match] => skip this observation
            for attr, name, typ in zip(attrs, attr_names, attr_types):
                if typ == 'eval':
                    val_1 = self.perform_eval(obs, attr)
                elif typ == 'name':
                    val_1 = self.extract_value(obs, attr)
                elif typ == 'pseudo':
                    val_1 = self.extract_pseudo(obs, nobs, attr)
                else:
                    assert False, 'Statement should not be reached.'
                vals[name].append(val_1)
        # return the dictionary-of-lists
        return vals

    def summarize(self, args, econo=True):
        r'''Find the summary of the sim as a dictionary held within the object.

        The convention is that the summary is built up by calling a series of analytic
        routines, each of which returns a dictionary of summary information.  The overall
        summary dictionary is a union of each individual summary.
        If econo, delete the DRM and keep only the summary.'''
        # this dict holds reductions for the current sim
        summary = self.extract_attrs(args.match, args.attrs, args.attr_names, args.attr_types)
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
    if verb > 1:
        print('Processing <%s> in pid #%d' % (f, os.getpid()))
    sim = SimulationRun(f)
    return sim.summarize(args)


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
            sims = list(map(SimulationRun, sim_files))
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
            # Note: for python3, need to ensure we materialize the list
            reductions = list(map_function(partial(outer_load_and_reduce, args=args, verb=self.args.verbose),
                                      self.sim_files))
        # hacky fix for no-drm case
        if len(self.sim_files) == 0:
            reductions = outer_load_and_reduce(None)
        # re-group the above reductions across sims
        # result is a dict containing reduced data, stored as self.summary
        self.regroup_and_accum(reductions, args.attr_names)

    def regroup_and_accum(self, reductions, attr_names):
        r'''Accumulate various summaries across the ensemble.

        Nothing is returned: result is placed in the object state.'''
        # flatten the reductions from [drm][attribute] to [attribute]
        # summary is a dictionary of lists
        summary = {}
        for attr in attr_names:
            summary[attr] = []
            for r in reductions:
                summary[attr].extend(r[attr])
        # record these summaries in the object
        self.summary = summary

    def dump(self, args, outfile):
        r'''Dump reduced data to output.'''
        saved_fields = args.attr_names[:]
        Nrow = len(self.summary[saved_fields[0]])
        # dump to CSV  or JSON
        if not args.json:
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
        else:
            # make an object, then dump it to JSON
            dumpable = []
            for i in range(Nrow):
                # dictionary mapping field -> value -- everything is a scalar here
                d = {f:self.summary[f][i] for f in saved_fields}
                dumpable.append(d)
            json.dump(dumpable, outfile, indent=2, cls=NumpyEncoder)


########################################
###
###  Mundane argument processing
###
########################################

def process_attr_program_inputs(args):
    r'''Process given attribute inputs, whether arguments or in a file.

    Inputs and outputs by the args.___ namespace.'''
    if args.file:
        # attributes-to-print are in a JSON file
        try:
            # ensure the load preserves order
            d = json.load(args.file, object_pairs_hook=collections.OrderedDict)
        except:
            sys.stdout.write('{}: Error loading json script from file "{}"\n'.format(args.progname, args.file.name))
            sys.stdout.write('{}: Traceback follows.\n'.format(args.progname))
            raise
        # close file; remember only the name so "args" can be serialized
        args.file.close()
        args.file = args.file.name
        assert isinstance(d, collections.abc.Mapping), 'JSON "{}" must translate to a dictionary'.format(args.file)
        # allow existing args.match to override
        if not args.match:
            args.match = d.pop('__match__', args.match)
        # update seed and obs_num from d if possible
        args.seed    = bool(d.pop('__seed__',    args.seed))
        args.obs_num = bool(d.pop('__obs_num__', args.obs_num))
        # update args.header if present in d
        args.header = bool(d.pop('__header__', args.header))
        # get a __name__ from d if possible, could be used for output
        args.name = d.pop('__name__', 'drm-attrs')
        # compose a list of strings to parse down into names and attributes
        #   ...as if it had been supplied by command line options
        #   ...omit _xxx attributes: private like comments, authorship, etc.
        attrs_given = ['{}:{}'.format(key, val) for key, val in d.items() if not key.startswith('_')]
    else:
        # attributes-to-print are in a list given on the command line
        attrs_given = args.attrs_orig
    # process the list of attributes to extract
    process_attr_list(args, attrs_given)
            

def process_attr_list(args, attrs_given):
    r'''Massage the list of attr arguments into matching lists
    that we will use later. Makes new fields in args using
    the attrs_given input dictionary, as follows:
      args.attrs        -- attribute specification
      args.attr_names   -- name of output attribute
      args.attr_types   -- type of attribute: 'name' or 'eval'
    '''
    def is_named(attr):
        r'''Return true if attr is an extended attribute name.
        This could be a one-liner RE, but left it this way for clarity.'''
        fields = attr.split('.')
        for f in fields:
            if re.match('-?[0-9]+$', f):
                continue # a numeric list index is named
            elif re.match('[a-zA-Z][a-zA-Z0-9_]*$', f):
                continue # all current drm entries are variable names
            # this field is not a dict-style drm entry, or a list index
            return False
        # all fields matched
        return True
        
    # start generated expressions (if any) at 1
    expr_num   = 1
    # allow comma-separated attr lists
    attrs      = []
    attr_names = []
    attr_types = []
    # handle pseudo-attributes
    if args.seed:
        attrs.append('seed')
        attr_names.append('seed')
        attr_types.append('pseudo')
    if args.obs_num:
        attrs.append('obs_num')
        attr_names.append('obs_num')
        attr_types.append('pseudo')
    # handle explicit attributes
    for attr_orig in attrs_given:
        # pre-scan to allow for comma-separated 'name' attributes
        # (issue here is that , might be in an 'eval' attribute)
        do_split = False
        if ',' in attr_orig and ':' not in attr_orig:
            # consider split if the attr_orig looks promising
            do_split = all(is_named(attr.split(',')))
        # split attributes into list if needed
        if do_split:
            attr_list = attr_orig.split(',')
        else:
            attr_list = [attr_orig]
        # iterate over each such attribute
        for attr in attr_list:
            # can happen if -s AND -a seed are given
            if attr in attrs:
                continue
            # pull away the field name, if it was given (name:...)
            prefix_found = re.match('[a-zA-Z][a-zA-Z0-9_]*\s*:', attr)
            if prefix_found:
                # save the name, and the actual attr
                attr_name = attr[:prefix_found.end()-1]
                attr = attr[prefix_found.end():]
            else:
                attr_name = ''
            # now we can decide if attr is name or eval type
            typ = 'name' if is_named(attr) else 'eval'
            # find name if needed
            if not attr_name:
                if typ == 'name':
                    attr_name = attr
                elif typ == 'eval':
                    attr_name = 'Expr_{}'.format(expr_num)
                    expr_num += 1
            # append all this information
            attrs.append(attr)
            attr_types.append(typ)
            attr_names.append(attr_name)
    # update args according to above
    args.attrs = attrs
    args.attr_names = attr_names
    args.attr_types = attr_types
    return


########################################
###
###  main()
###
########################################

def main(args, outfile):
    if VERBOSITY:
        sys.stderr.write('%s: Loading %d file patterns.\n' % (args.progname, len(args.infile)))
    lazy = True
    ensemble = EnsembleSummary(args.infile, args, lazy=lazy)
    if ensemble.Ndrm_actual == 0:
        print('%s: Warning: no actual DRMs present.' % (args.progname, ))
    # save a reference to the universe (stars and their attributes)
    # don't want to assume a SPC
    # ensemble.spc = ensemble.sims[0].spc
    # ensemble.Nstar = ensemble.sims[0].Nstar
    if VERBOSITY:
        sys.stderr.write('%s: Tabulating.\n' % args.progname)
    ensemble.load_and_reduce()
    # print '%s: Dumping.' % args.progname
    # default is ~70, which means small-sized vectors are line-wrapped in the CSV
    np.set_printoptions(linewidth=1000)
    ensemble.dump(args, outfile)


########################################
###
###  Program Entry
###
########################################

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Extract attributes from DRMs and send to stdout.",
                                     epilog='')
    # either -a or -f must be given, not both
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-a', '--attr', type=str, action='append', dest='attrs_orig', metavar='ATTR',
                           help='Repeat -a for more attributes')
    group.add_argument('-f', '--file', type=argparse.FileType('r'), help='JSON file of attributes')
    parser.add_argument('drm', metavar='DRM', nargs='+', default=[], help='drm file list')
    parser.add_argument('-m', '--match', type=str, default='',
                            help='Attribute to match')
    parser.add_argument('-v', help='Verbosity', action='count', dest='verbose',
                            default=0)
    parser.add_argument('-1', help='Supply header line', action='store_true', dest='header',
                            default=False)
    parser.add_argument('--json', help='JSON output format', action='store_true', dest='json',
                            default=False)
    parser.add_argument('-d', '--delimiter', help='CSV field delimiter', type=str)
    parser.add_argument('-s', '--seed', help='Output the DRM seed', action='store_true')
    parser.add_argument('-n', '--obs_num', help='Output the observation number', action='store_true')
    parser.add_argument('-j', '--jobs', help='Number of parallel jobs (default = %(default)d)',
                      type=int, dest='jobs', default=int(N_CPU*0.65))
    args = parser.parse_args()
    
    # set umask in hopes that files will be group-writable
    os.umask(0o002)

    VERBOSITY = args.verbose

    # program name, for convenience
    args.progname = os.path.basename(sys.argv[0])
    
    # do it like this for now
    args.infile = args.drm

    # process attr args here, before calling main
    #   uses: args.attrs_orig or attrs.file [which attrs]
    #   sets: args.attrs
    #         args.attr_names
    #         args.attr_type
    process_attr_program_inputs(args)
    if VERBOSITY > 0:
        sys.stderr.write('Attributes to be retrieved:\n')
        for attr, name, typ in zip(args.attrs, args.attr_names, args.attr_types):
            sys.stderr.write('  retrieving "{}" type attribute: {} <- {}\n'.format(typ, name, attr))

    # no reason to complicate this for now
    outfile = sys.stdout

    main(args, outfile)
    sys.exit(0)
