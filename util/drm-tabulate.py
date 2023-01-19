#!/usr/bin/env python
r'''
drm-tabulate.py: extract named fields from a pile of DRMs

usage:
  drm-tabulate.py [-1sn] [--json] [ -m ATTR ] [-a ATTR | -f FILE] DRM ...

where:
  DRM ... is a list of DRM pickles

Choose attributes of an observation "obs" in the DRM with:
  -a ATTR -> output obs[ATTR] (repeat -a OK, see below)
  -S ATTR -> output the SPC attribute ATTR for obs[star_ind]
  -P ATTR -> output the SPC attributes ATTR for obs[plan_inds[:]]
  -p ATTR -> output the list of SPC attributes ATTR for obs[plan_inds]
  -e CODE -> compute a value by [e]valuating the Python CODE expression

Match a subset of observations with:
  -m ATTR -> produce output if ATTR is a keyword in the observation

Output special pseudo-attributes:
  -s: output the DRM seed number (__seed__, boolean)
  -n: output the DRM observation number (__obs_num__, boolean)
  --plan_num: output the planet number (__plan_num__, boolean)
  -1: print line 1, the CSV header (__header__, boolean)

Take the above options from a file with:
  -f FILE: take ATTRs from FILE *instead of* -a/-S/-P/-p options (see below)

Some less-useful options:
  -A: list available DRM and SPC attributes and values on stderr, for reference
  --json: output is JSON, rather than standard CSV
  --empty: output a record when planet attributes selected, even if no planets present
  -v: increase verbosity (output is to stderr)
  -j N: use N parallel workers (default = ~2/3 of cores)
        If N = 0 or 1, no parallelism: needed for debugging.

DRM ATTRIBUTES
--------------

The main degree of freedom is attribute specification -- which is done
using either of two notations. Below, suppose "obs" is one entry in the DRM.

-a => Directly named attributes
    1) Get a field in obs by giving: "-a arrival_time"
    2) Drill into nested attributes with "."; for example
          obs['char_mode']['lam']
        is extracted with:
          -a char_mode.lam
        For lists, give the index number; obs['plan_inds'][0] is
          -a plan_inds.0 
        The last phase angle, obs['char_params'][-1], is:
          -a char_params.phi.-1
    3) Supply a comma-separated list of such DRM fields in one -a
          -a arrival_time,slew_time,scMass

-e => Evaluated expressions
    1) A Python expression can be given, which is evaluated in the context 
        of variables named for each field in "obs". To output a count of 
        detections using the obs['det_status'] list, use:
          -e "np.sum(det_status == 1)"
    2) No comma-separated lists are allowed, due to ambiguity.
    3) The full SPC is available, if desired, using spc[...], so
       -e "spc['Spec'][star_ind]"  <==> -S Spec
       See below for more on -S.

For either -a or -e, the resulting column can be custom-named with 
a label:attr construct, such as
          -a "lambda:char_mode.lam"
          -a "det_count:np.sum(det_status == 1)"
otherwise a basic generated name is used.
(But: Attributes in comma-separated expressions cannot be custom-named.)

STAR-PLANET ATTRIBUTES
----------------------

-S => shortcut for Star attributes
    -S ATTR means: look up the named ATTR for obs['star_ind'] in
    the corresponding SPC file, e.g.
      -S Spec => spc['Spec'][obs['star_ind']]
    means to output the spectral class of obs['star_ind']

-P => shortcut for Planet attributes
    -P ATTR means to look up the named ATTR for each planet in 
    obs['plan_inds'] in the SPC file, e.g.
      -P Mp => spc['Mp'][obs['plan_inds']]
    Note that obs['plan_inds'] is in general a vector. 
    So, if you give -P, this program "scalar-expands" the vector
    to write one row of output *for each plan_ind in plan_inds*. 
    This facilitates row-by-row processing.
    Note that if plan_inds = [], no record will be written. 
    To write a record in the zero-planet case anyway, specify --empty.

-p => shortcut for alternate Planet attributes
    This is the same lookup as -P, but the vector is output to
    that *one* column in the row. This would be more useful for 
    the JSON output; the CSV format looks like:
       "[0.282, 1.044]"

For all of the above, an optional column-name can be given just as
for -a and -e. If not given, the attribute name, or a generated
string, will be used.

The SPC file is loaded using the filename convention that 
  .../drm/NAME.pkl -> .../spc/NAME.spc
If no -S/-P/-p is given, the SPC is not loaded, to allow use of 
this program when only the DRM is present. If the SPC file is needed
to support spc[...] within "-e" constructs above, load of the SPC 
can be forced by giving -S "".


EXTERNAL FILE
-------------

These expressions can become complex, so the "-a ATTR" and all above
constructs can be placed in a JSON file and specified with -f FILE:

{
  "__match__": "char_info",
  "lambda": "char_mode.lam",
  "__eval__": {
    "char_count": "np.sum(char_info[0]['char_status'] == 1)"
  },
  "__star__": {
    "Spec": "Spec",
    "Name": "Name"
  }
}

Recall JSON uses double-quotes for strings, and Python dictionary keys
are extracted with ...['attr'], so there is no conflict between quote marks.
The key name (e.g., "char_count" above) is used as the column name for the CSV.
Above, "__match__" abbreviates the -m construct, so the JSON file can be
self-contained. The full list is:
   "__eval__"    -> -e 
   "__star__"    -> -S
   "__planet__"  -> -P  (one row per planet)
   "__planets__" -> -p  (list placed in one field)
Other program flags (like -s) can be given by Booleans in the 
JSON file as well. See the top of this usage note for the __attribute_name__
controlling each flag.


USAGE
-----

Typical usage:
  # each arrival_time
  util/drm-tabulate.py -a arrival_time sims/HabEx_4m_dmag26/drm/*

  # each char_time, looking only at chars
  util/drm-tabulate.py -m char_time -a char_time sims/HabEx_4m_dmag26/drm/*

  # number of successful chars, only for chars
  util/drm-tabulate.py -m char_info -a "np.sum(char_info[0]['char_status'] == 1)" sims/HabEx_4m_dmag26/drm/*

  # planet-by-planet output of: char_status, planet mass, star spectral class, etc.
  util/drm-tabulate.py -ns1 -m char_time -e "CS:char_status" -P Mp -e "SpecLetter:[spc['Spec'][star_ind][0]]" -S Spec -a ct:char_time -a char_mode.lam sims/.../drm/*.pkl

turmon apr 2019, feb 2022, dec 2022
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
import pprint
import warnings
import collections
from functools import partial
from dataclasses import dataclass
from typing import ClassVar
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

# used to limit multiple attribute-printing under -A
RECORD_WAS_SHOWN = False

########################################
###
###  Utility Functions
###
########################################

# valid attribute flavors (for input args)
AttrFlavors = ['pseudo', 'attr', 'eval', 'star', 'planet', 'planets']

# container for a single Attribute (field to tabulate)
@dataclass
class Attribute:
    names_used: ClassVar = set()
    flavor: str
    name: str
    text: str
    # raises if name is duplicated -- better now than later
    def __post_init__(self):
        assert self.flavor in set(AttrFlavors), 'Unrecognized flavor <{}>'.format(self.flavor)
        assert self.name not in self.names_used, 'Duplicate name <{}>'.format(self.name)
        self.names_used.add(self.name)


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
        # load a spc file - only if needed to lookup -S, -P attrs
        if args.load_spc:
            g = f.replace('pkl', 'spc').replace('/drm/', '/spc/')
            if os.path.isfile(g):
                spc = pickle.load(open(g, 'rb'), **PICKLE_ARGS)
            else:
                raise ValueError('Could not find a .spc file to match DRM <%s>' % f)
        else:
            spc = None
        gc.enable()
        # set up object state
        self.name = f
        self.seed = int(os.path.splitext(os.path.basename(f))[0])
        # self.Nstar = len(spc['Name'])
        self.spc = spc # = None, if SPC not needed
        self.drm = drm
        self.summary = None # place-holder
    

    def extract_value_obs(self, obs, attr, just_peeking=False):
        r'''Extract a value, attr, from an observation, obs, allowing recursive lookup.'''
        if not attr:
            # recursion base case:
            #   our obs = caller's obs[seg0] was the last lookup,
            #   and our attr = caller's segments[1:] is empty
            if just_peeking:
                # original obs[attr] is indeed in obs
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
                pp = pprint.PrettyPrinter(indent=4)
                # announce the cause now, exception will raise below
                print('Fatal: Value extraction failed during {}.'.format(
                    'list indexing' if type(seg0) is int else 'dict lookup'), file=sys.stderr)
                print('No "{}" found in obs object (or sub-object). Attribute naming error?'.format(seg0), file=sys.stderr)
                print('Object or sub-object ({}) = '.format(type(obs)), file=sys.stderr)
                print(pp.pformat(obs), file=sys.stderr)
            # recursive extraction from obs[seg0]
            sub_object = obs[seg0]
            return self.extract_value_obs(sub_object, '.'.join(segments[1:]), just_peeking=just_peeking)

    def extract_pseudo(self, obs, attr, nobs=None):
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
        elif attr == 'plan_num':
            return 0 # likely altered later, depending on #planets
        else:
            assert False, 'Should not be reached (attr = {})'.format(attr)
            
    def perform_eval(self, obs, attr):
        r'''Extract attribute by evaluating attr in the context of obs.'''
        try:
            val = eval(attr, {"np": np, "spc": self.spc, **obs})
        except:
            sys.stderr.write('Error: Failed to eval() the given attribute.\n')
            sys.stderr.write('Evaluated string is surrounded by || below:\n|{}|\n'.format(attr))
            sys.stderr.write('Context is:\n')
            sys.stderr.write(repr(obs) + '\n')
            raise
        return val

    def extract_value_star(self, obs, attr):
        r'''Extract star-related attribute by evaluating attr in the context of obs.'''
        assert self.spc, 'No SPC info loaded for lookup of {}\n'.format(attr)
        try:
            sind = obs['star_ind']
            val = strip_units(self.spc[attr][sind])
        except:
            print('Error: Did not find "{}" in SPC for star index {}.'.format(attr, sind), file=sys.stderr)
            print('Properties available for -S include:', file=sys.stderr)
            print('\t' + '\n\t'.join(self.spc.keys()), file=sys.stderr)
            print('Consider using -A for more information.', file=sys.stderr)
            raise
        return val

    def extract_value_planet(self, obs, attr):
        r'''Extract planet-related attribute by evaluating attr in the context of obs.'''
        assert self.spc, 'No SPC info loaded for lookup of {}\n'.format(attr)
        try:
            plan_inds = obs['plan_inds']
            val = [strip_units(self.spc[attr][p]) for p in plan_inds]
        except:
            print('Error: Did not find "{}" in SPC for planet index {}.'.format(attr, repr(plan_inds)), file=sys.stderr)
            print('Properties available for -P/-p include:', file=sys.stderr)
            print('\t' + '\n\t'.join(self.spc.keys()), file=sys.stderr)
            print('Consider using -A for more information.', file=sys.stderr)
            raise
        # val is always a list from this function, we'll process further later
        return val

    def extract_value_any(self, obs, attr, **kwargs):
        r'''Selector function that switches on the attribute flavor to extract a value.'''
        flavor, text = attr.flavor, attr.text
        if flavor == 'attr':
            val = self.extract_value_obs(obs, text)
        elif flavor == 'eval':
            val = self.perform_eval(obs, text)
        elif flavor == 'star':
            val = self.extract_value_star(obs, text)
        elif flavor == 'planet':
            # val is a length=#planets list
            val = self.extract_value_planet(obs, text)
        elif flavor == 'planets':
            # value (which is expected to already be a list here)
            # is boxed into a len = 1 list
            val = self.extract_value_planet(obs, text)
        elif flavor == 'pseudo':
            val = self.extract_pseudo(obs, text, **kwargs)
        else:
            assert False, 'Statement should not be reached: unrecognized flavor'
        # convert np.ndarrays -> lists (so we can use .extend later on any returned value)
        # (but: ndim == 0 -> scalar -> leave it)
        if isinstance(val, np.ndarray) and val.ndim > 0:
            val = val.tolist()
        # box everything but 'planet' and 'eval' into a len = 1 list
        # this len=1 list may be extended later. Here are the exceptions:
        #  planet -> returns a list anyway (b/c spc[attr][plan_inds] is a list)
        #  eval -> if we /always/ box "val" here, a list returned by eval will be 
        #          double-boxed and can't be extended later. OTOH, if boxing 
        #          a non-planet result is desired, eval can be forced to return 
        #          a list by surrounding in [...]
        if flavor not in ('planet', 'eval'):
            val = [val]
        # - the code below recognizes Py/NP vectors and separates them from
        # strings (which also have a len), so we can vectorize correctly later
        # - it prevents returning non-boxed values, because values are vectorized later
        # with val.extend(...)
        # - correct behavior is that -P Mp and -e char_status both vectorize, but
        # -S Spec (string) does not.
        if True:
            # to vectorize later, we need the length. length=1 if scalar or string.
            # else, allow for lists that are np vectors or python lists
            try:
                # this can be a boxed list of length=1
                val_len = len(val)
            except TypeError:
                # believe this triggers only on scalars from -e, must always box val once
                val = [val]
                val_len = 1
            if isinstance(val, str):
                val_len = 1 # str alone is OK for py3 or np.str
        return val_len, val

    def show_attributes(self, obs):
        r'''Pretty-print available attributes to stderr.'''
        # global-var synchronization scheme only works with -j 1: adequate
        global RECORD_WAS_SHOWN
        if RECORD_WAS_SHOWN:
            return
        RECORD_WAS_SHOWN = True
        pp = pprint.PrettyPrinter(indent=4)
        # obs attributes
        print('DRM observation attributes [use -a]:', file=sys.stderr)
        # loop over keys => attributes appear as you'd name them with -a
        for k in sorted(obs.keys()):
            value_pp = pp.pformat(obs[k])
            print('  {}: {}'.format(k, value_pp), file=sys.stderr)
        if self.spc:
            print('', file=sys.stderr)
            print('SPC attributes [use -S or -P]:', file=sys.stderr)
            # loop over keys => attributes appear as you'd name with -S/-P
            for k in sorted(self.spc.keys()):
                try:
                    xtra = 'array' + str(self.spc[k].shape)
                    if self.spc[k].size < 10:
                        # if it's short, give its values
                        xtra = xtra + ' = ' + pp.pformat(self.spc[k])
                    else:
                        # otherwise, give its first few values only
                        num_show = min(self.spc[k].size, 3)
                        xtra = xtra + ' = ' + ', '.join([pp.pformat(x) for x in (self.spc[k][:num_show])])
                        if num_show < self.spc[k].size:
                            xtra = xtra + ', ...'
                except AttributeError:
                    # (no .shape -> not numpy)
                    xtra = pp.pformat(self.spc[k])
                print('  {}: {}'.format(k, xtra), file=sys.stderr)
        else:
            print('Note: SPC not loaded. Use -S "" to force load.', file=sys.stderr)


    def extract_attrs(self, match, attrs, show_attrs):
        r'''Extract attributes from each obs in the DRM.

        Returns a dictionary mapping attributes to lists, one list entry
        per DRM observation.
        '''
        # accumulate values in these lists, one per attribute
        vals = {name: [] for name in attrs.keys()}
        not_yet_shown = True
        for nobs, obs in enumerate(self.drm):
            # 1: no obs[match] => skip this observation
            if match and not self.extract_value_obs(obs, match, just_peeking=True):
                continue
            if show_attrs and not_yet_shown:
                not_yet_shown = False
                self.show_attributes(obs)
            # 2: extract all needed values -- note, vals_1 is set for every name
            val_lens, vals_1 = dict(), dict()
            for name, attr in attrs.items():
                val_lens[name], vals_1[name] = self.extract_value_any(obs, attr, nobs=nobs)
            # 3A: detect #planets
            planet_guesses = set(val_lens.values()) - set([1])
            assert len(planet_guesses) <= 1, 'Two attributes have incommensurate non-scalar lengths: should not happen'
            if planet_guesses:
                Nplan = planet_guesses.pop()
            else:
                Nplan = 1
            ### print('EX_ATTR: nobs = {}, Nplan = {}'.format(nobs, Nplan))
            # 3B: extend vals_1 along planets if needed
            # 3B.1 -- Nplan = 0 case
            #   skip this record if there was a planet attribute selected, but no planets there
            #   otherwise, we continue and the PlanetNum will tabulate as 0
            if Nplan == 0:
                if args.empty_skipped:
                    continue
                else:
                    # Nplan == 0 case, but not empty_skipped
                    # ==> we will finish the loop, and output a record for this obs
                    # fields will contain whatever extract_value() returned above
                    # plan_num, if requested, will be 0
                    pass
            # 3B.2 -- Nplan >= 1 case
            #   pad the other fields downward to match the planet vector
            #   set plan_num field if needed -- [1:Nplan], excluding 0
            #   test is ">= 1" here to set plan_num for Nplan == 1
            if Nplan >= 1:
                for name in vals_1.keys():
                    vals_1[name].extend([vals_1[name][-1]] * (Nplan - val_lens[name]))
                if 'plan_num' in vals_1:
                    # set plan_num correctly ... its bogey value is [0]
                    vals_1['plan_num'] = [(index + 1) for index in range(Nplan)]
            # 4: vals += vals_1 (both are lists)
            for name in vals_1.keys():
                vals[name].extend(vals_1[name])
        # return the dictionary-of-lists
        return vals

    def summarize(self, args, econo=True):
        r'''Find the summary of the sim as a dictionary held within the object.

        The convention is that the summary is built up by calling a series of analytic
        routines, each of which returns a dictionary of summary information.  The overall
        summary dictionary is a union of each individual summary.
        If econo, delete the DRM and keep only the summary.'''
        # this dict holds reductions for the current sim
        summary = self.extract_attrs(args.match, args.all_attrs, args.show_attributes)
        # delete the base data if asked
        if econo:
            self.drm = None
            self.spc = None
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
    def __init__(self, in_files, args):
        '''Load an ensemble of simulations.'''
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
        self.regroup_and_accum(reductions, args.all_attrs.keys())

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
        saved_fields = args.all_attrs.keys()
        Nrows = set(len(self.summary[fname]) for fname in saved_fields)
        if len(Nrows) != 1:
            warning('Attribute lists have differing lengths') # shouldn't happen
        Nrow = max(Nrows)
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

def json_to_object(args):
    r'''Load a JSON dictionary, extract any program flags, return the remainder.'''
    # 1: attempt to load a dict from a JSON file
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
    # 2: Propagate simple flags from d -> args
    # allow existing args.match to override
    if not args.match:
        args.match = d.get('__match__', args.match)
    # update args.header if present in d
    args.header = bool(d.get('__header__', args.header))
    # 3: Make containers for __attr__ and __pseudo__
    # Insert container duplicating top-level attributes, for later
    d['__attr__'] = {key:val for key, val in d.items() if not key.startswith('_')}
    # Insert container of each pseudo-attribute that is present and truthy
    d['__pseudo__'] = {key:key for key in ('seed', 'obs_num', 'plan_num') if d.get('__{}__'.format(key))}
    return d


def argtexts_to_object(args):
    r'''Convert lists of argument-texts to a dict-of-dicts indexed by flavor.'''

    # make a unique name for flavor (Eval_1, Attr_5, etc.)
    count = collections.Counter() # persistent state for below
    def generated_name(flavor, count=count):
        count[flavor] += 1
        return '{}_{}'.format(flavor.capitalize(), count[flavor])

    # iterate over attribute-lists from all sources (-a, -e, etc.)
    obj = {}
    for flavor in AttrFlavors:
        obj[flavor] = {}
        # 1: handle comma-splits
        # iterate over each supplied textual attribute
        # (-a foo,bar -a baz => argtext='foo,bar', argtext='baz')
        argtexts = []
        for argtext in vars(args)[flavor]:
            # split argtext on commas, unless it was eval (-e)
            if flavor != 'eval':
                argtext_split = argtext.split(',')
            else:
                argtext_split = [argtext]
            argtexts.extend(argtext_split)
        # 2: extract name: prefix if present; if not, make a name
        for argtext in argtexts:
            # pull away the field name, if it was given (name:...)
            # pattern: whitespace(name)whitespace:whitespace(TEXT)whitespace
            name_colon_value = re.search(r'\s*(\w+)\s*:\s*(.*\S)\s*', argtext)
            if name_colon_value:
                # save the name, and the actual attr
                attr_name, attr_text = name_colon_value.groups()
            elif flavor in ('pseudo', 'attr', 'star', 'planet'):
                attr_name, attr_text = argtext, argtext
            elif flavor == 'planets':
                # otherwise, attribute name-clash (-P/-p) is invited
                attr_name, attr_text = argtext + '_list', argtext
            else:
                # will generate a name below
                attr_name, attr_text = generated_name(flavor), argtext
            obj[flavor][attr_name] = attr_text
    return obj


def process_attr_program_inputs(args):
    r'''Process given attribute inputs, whether arguments or in a file.

    The args.___ namespace is updated by the JSON file contents.
    The return value is a dict of Attribute objects, indexed
    by the attribute name ("column name") in the output tablulation.'''
    # Strategy: under either the JSON-file or command-line setup, convert
    # attribute specification into a common-denominator object: a dict-of-dicts
    # mapping flavor -> {name1:text1, name2:text2, ...}
    if args.file:
        # dictionary of the JSON file
        d = json_to_object(args)
        # make dict-of-dicts object
        obj = {}
        for flavor in AttrFlavors:
            obj[flavor] = d.get('__{}__'.format(flavor), {})
    else:
        # attributes-to-print are in argument-lists from the command line
        obj = argtexts_to_object(args)
    # Compose the list of all attributes we will extract (across all
    # flavors) by collecting attributes over flavors and names
    # Output field order preserved b/c dict is ordered (py 3.7+)
    all_attrs = dict()
    for flavor, attr_dict in obj.items():
        for name, attr_text in attr_dict.items():
            all_attrs[name] = Attribute(flavor, name, attr_text)
    return all_attrs


########################################
###
###  main()
###
########################################

def main(args, outfile):
    if VERBOSITY:
        sys.stderr.write('%s: Loading %d file patterns.\n' % (args.progname, len(args.infile)))
    ensemble = EnsembleSummary(args.infile, args)
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
    group = parser.add_argument_group('Attributes via command line')
    group.add_argument('-a', '--attr', type=str, action='append', default=[], metavar='ATTR',
                           dest='attr', help='Named attribute within obs (e.g., arrival_time)')
    group.add_argument('-e', '--eval', type=str, action='append', default=[], metavar='CODE',
                           dest='eval', help='Evaluate expression within obs namespace')
    group.add_argument('-S', '--star', type=str, action='append', default=[], metavar='ATTR',
                            dest='star', help='Star attribute name within SPC via sind, e.g. Spec')
    group.add_argument('-P', '--planet', type=str, action='append', default=[], metavar='ATTR',
                            dest='planet', help='Planet attribute within SPC via plan_inds, e.g. Mp')
    group.add_argument('-p', '--planets', type=str, action='append', default=[], metavar='ATTR',
                            dest='planets', help='Planet attribute, as above, as a list per obs')
    group = parser.add_argument_group('Pseudo-attributes')
    # -s and -n put a pseudo-attribute into that list
    group.add_argument('-s', '--seed', dest='pseudo', help='Output the DRM seed', action='append_const', const='seed', default=[])
    group.add_argument('-n', '--obs_num', dest='pseudo', help='Output the observation number', action='append_const', const='obs_num')
    group.add_argument('--plan_num', dest='pseudo', help='Output the planet number', action='append_const', const='plan_num')
    group = parser.add_argument_group('Attributes by JSON file')
    group.add_argument('-f', '--file', type=argparse.FileType('r'), help='JSON file of attributes')

    parser.add_argument('drm', metavar='DRM', nargs='+', default=[], help='drm file list')
    parser.add_argument('-m', '--match', type=str, default='', metavar='ATTR',
                            help='Attribute within obs to match, else, skip the obs')
    group = parser.add_argument_group('Output format')
    group.add_argument('-1', help='Supply CSV header line', action='store_true', dest='header',
                            default=False)
    group.add_argument('--delimiter', help='CSV field delimiter', metavar='STR', type=str)
    group.add_argument('--json', help='JSON output format', action='store_true', dest='json',
                            default=False)
    group = parser.add_argument_group('Seldom used')
    group.add_argument('--empty', help='output records despite having 0 planets', action='store_false', dest='empty_skipped',
                            default=True)
    group.add_argument('-A', help='Show list of Available Attributes on stderr', action='store_true', dest='show_attributes')
    group.add_argument('-v', help='Verbosity', action='count', dest='verbose', default=0)
    group.add_argument('-j', '--jobs', help='Number of parallel jobs (default = %(default)d)',
                      type=int, dest='jobs', default=int(N_CPU*0.65))
    args = parser.parse_args()
    
    # set umask in hopes that files will be group-writable
    os.umask(0o002)

    VERBOSITY = args.verbose

    # program name, for convenience
    args.progname = os.path.basename(sys.argv[0])
    
    # do it like this for now
    args.infile = args.drm

    # do we need to load the SPC file?
    args.load_spc = args.star or args.planet or args.planets

    # enforce that -A => -j 1
    if args.show_attributes:
        args.jobs = 1

    # Enforce that -P => --plan_num
    if args.planet and 'plan_num' not in args.pseudo:
        args.pseudo.append('plan_num')

    # process attr args here, before calling main
    args.all_attrs = process_attr_program_inputs(args)
    if VERBOSITY > 0:
        sys.stderr.write('Attributes to be retrieved:\n')
        for name, attr in args.all_attrs.items():
            sys.stderr.write('  retrieving "{}" type attribute: {} <- {}\n'.format(attr.flavor, name, attr.text))

    # no reason to complicate this for now
    outfile = sys.stdout

    main(args, outfile)
    sys.exit(0)
