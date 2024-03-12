#!/usr/bin/env python
r'''
drm-tabulate.py: extract named fields from a pile of DRMs

usage:
  `drm-tabulate.py [-1sn] [--json] [ -m ATTR ] [-a ATTR | -f FILE] DRM ...`

Args:
  DRM (file): a list of DRM pickles

Attributes of an observation "obs" in the DRM are printed by 
naming them in one of the following ways.

Attribute naming:
  + `-a ATTR`: output `obs[ATTR]` (repeat -a OK, see below)
  + `-S ATTR`: output the SPC attribute ATTR for `obs[star_ind]`
  + `-P ATTR`: output the SPC attributes ATTR for `obs[plan_inds[:]]`
  + `-p ATTR`: output the list of SPC attributes ATTR for `obs[plan_inds]`
  + `-e CODE`: compute a value by [e]valuating the Python CODE expression

Pseudo-attributes:
  + `-s`: output the DRM seed number (`__seed__`, boolean)
  + `-n`: output the DRM observation number (`__obs_num__`, boolean)
  + `--plan_num`: output the planet number (`__plan_num__`, boolean)
  + `-1`: supply the CSV header on line 1 (`__header__`, boolean)

Match only a subset of observations:
  + `-m ATTR`: produce output only if ATTR is a keyword in obs
  + `-M ATTR`: produce output if ATTR is NOT a keyword in obs

External file:
  + -f FILE: take ATTRs from FILE *instead of* `-a/-S/-P/-p` options (see below)

Less-useful options:
  + `-A`: list available DRM and SPC attributes and values on stderr, as a reference.
      It honors match (`-m`), and inverse match (`-M`); supply `-S ""` to get SPC attributes.
  + `--json`: output is JSON, rather than standard CSV
  + `--format`: supply a printf-style format string (e.g., %.3g) for CSV floats
  + `--empty`: output a record when planet attributes selected, even if no planets present
  + `-v`: increase verbosity (output is to stderr)
  + `-j N`: use N parallel workers (default = ~2/3 of cores)
        If N = 0 or 1, no parallelism: needed for debugging.

DRM ATTRIBUTES
--------------

The main degree of freedom is attribute specification -- which is done
using either of two notations. Below, suppose "obs" is one entry in the DRM.

* `-a` => Directly named attributes
    1. Get a field in obs by giving: `-a arrival_time`
    2. Drill into nested attributes with `.`. For example
          `obs['char_mode']['lam']`
        is extracted with:
          `-a char_mode.lam`
        For lists, give the index number; `obs['plan_inds'][0]` is
          `-a plan_inds.0`
        The last phase angle, `obs['char_params'][-1]`, is:
          `-a char_params.phi.-1`
    3. Supply a comma-separated list of such DRM fields at once with
          `-a arrival_time,slew_time,scMass`

* `-e` => Evaluated expressions
    1. A Python expression can be given, which is evaluated in the context 
        of variables named for each field in "obs". To output a count of 
        detections using the `obs['det_status']` list, use:
          `-e "np.sum(det_status == 1)"`
    2. No comma-separated lists are allowed, due to ambiguity.
    3. The full SPC is available, if desired, using `spc[...]`, so
       `-e "spc['Spec'][star_ind]`  <==> `-S Spec`
       See below for more on `-S`.

For either `-a` or `-e`, the resulting column can be custom-named with 
a `label:attr` construct, such as
``` shell
          -a "lambda:char_mode.lam"
          -e "det_count:np.sum(det_status == 1)"
```
otherwise a basic generated name is used.
(But: Attributes in comma-separated expressions cannot be custom-named.)

STAR-PLANET ATTRIBUTES
----------------------

+ `-S` => shortcut for Star attributes.

    `-S ATTR` means: look up the named ATTR for `obs['star_ind']` in
    the corresponding SPC file, e.g.
      `-S Spec` => `spc['Spec'][obs['star_ind']]`
    which will output the spectral class of `obs['star_ind']`

+ `-P` => shortcut for Planet attributes.

    `-P ATTR` means to look up the named ATTR for each planet in 
    `obs['plan_inds']` in the SPC file, e.g.
      `-P Mp` => `spc['Mp'][obs['plan_inds']]`.
    
    Note that `obs['plan_inds']` is in general a vector. 
    So, if you give `-P`, this program "scalar-expands" the vector
    to write one row of output *for each plan_ind in plan_inds*. 
    This facilitates row-by-row processing.

    Note that if `plan_inds = []`, no record will be written. 
    To write a record in the zero-planet case anyway, specify `--empty`.

+ `-p` => shortcut for alternate Planet attributes.

    This is the same lookup as `-P`, but the vector is output to
    that *one* column in the row. This would be more useful for 
    the JSON output; the CSV format looks like:
       `"[0.282, 1.044]"`
    which would not look like a number to downstream consumers.

For all of the above, an optional column-name can be given just as
for `-a` and `-e`. If not given, the attribute name, or a generated
string, will be used.

The SPC file is loaded using the filename convention that 
  `.../drm/NAME.pkl` goes with `.../spc/NAME.spc`
If no `-S/-P/-p` is given, the SPC is not loaded, to allow use of 
this program when only the DRM is present. If the SPC file is needed
to support `spc[...]` within `-e` constructs above, load of the SPC 
can be forced by giving `--load_spc` (or if using External File,
`"load_spc": true`).


EXTERNAL FILE
-------------

These expressions can become complex, so the "-a ATTR" and all above
constructs can be placed in a JSON file and specified with `-f FILE`:

``` json
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
```

Recall JSON uses double-quotes for strings, and Python dictionary keys
can be extracted with `d['attr']`, so there is no conflict between quote marks.
The column name (e.g., `char_count` above) is used as the column name for the CSV.
Above, `__match__` abbreviates the `-m` construct, so the JSON file can be
self-contained. The full list is:
```
   "__eval__"    -> -e 
   "__star__"    -> -S
   "__planet__"  -> -P  (one row per planet)
   "__planets__" -> -p  (list placed in one field)
   "__match_inv__"  -M  (inverse match)
```
Other program flags (like `-s`) can be given by Booleans in the 
JSON file as well. See the top of this usage note for the `__attribute_name__`
controlling each flag.

As on the command line, the SPC file is loaded if `__star__`, `__planet__`, 
or `__planets__ `is present. To force the load if an `__eval__` construct needs
the SPC file, use `--load_spc`, or the Boolean directive `"__load_spc__": true`
within the JSON.


USAGE
-----

Typical usage:
``` shell
  # each arrival_time
  util/drm-tabulate.py -a arrival_time sims/HabEx_4m_dmag26/drm/*

  # each char_time, looking only at chars
  util/drm-tabulate.py -m char_time -a char_time sims/HabEx_4m_dmag26/drm/*

  # number of successful chars, only for chars
  util/drm-tabulate.py -m char_info -e "np.sum(char_info[0]['char_status'] == 1)" sims/HabEx_4m_dmag26/drm/*

  # planet-by-planet output of: char_status, planet mass, star spectral class, etc.
  util/drm-tabulate.py -ns1 -m char_time -e "CS:char_status" -P Mp -e "SpecLetter:[spc['Spec'][star_ind][0]]" -S Spec -a ct:char_time -a char_mode.lam sims/.../drm/*.pkl
```

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
    flavor: str
    name: str
    text: str
    expand: bool
    def __post_init__(self):
        assert self.flavor in set(AttrFlavors), f'Unrecognized flavor <{self.flavor}>'
        assert '+' not in self.name, f'Illegal character in attribute name {self.name}'


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
    def __init__(self, f, load_spc=False):
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
        if load_spc:
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
        r'''Extract a value, attr, from an observation, obs, allowing recursive lookup.

        The just_peeking argument is to detect presence of an attribute.
        It returns True/False, not the attribute itself.'''
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
            pp = pprint.PrettyPrinter(indent=4)
            print('Error: Failed to eval() the given attribute.', file=sys.stderr)
            print('String we attempted to eval() is below inside ||:', file=sys.stderr)
            print('\t|{}|'.format(attr), file=sys.stderr)
            self.show_attributes(obs, msg_obs='eval context is:', msg_spc=' spc: {')
            print(' }', file=sys.stderr)
            print('\nTraceback follows.', file=sys.stderr)
            # print(pp.pformat(context), file=sys.stderr)
            raise
        # for parallelism with other value-getters, strip the units now
        return strip_units(val)

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
            # record that we boxed it already
            is_list = True
        else:
            is_list = False
        # box everything but 'planet' and 'eval' into a len = 1 list
        # this len=1 list may be extended later. Here are the exceptions:
        #  planet -> returns a list anyway (b/c spc[attr][plan_inds] is a list)
        #  eval -> if we /always/ box "val" here, a list returned by eval will be 
        #          double-boxed and can't be extended later. OTOH, if boxing 
        #          a non-planet result is desired, eval can be forced to return 
        #          a list by surrounding in [...]
        if flavor not in ('planet', 'eval') and not is_list:
            val = [val]
        # - the code below recognizes Py/NP vectors and separates them from
        # strings (which also have a len), so we can vectorize correctly later
        # - it prevents returning non-boxed values, because values are vectorized later
        # with val.extend(...)
        # - correct behavior is that -P Mp and -e char_status both vectorize, but
        # -S Spec (string) does not.

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

    def show_attributes(self, obs, msg_obs=None, msg_spc=None):
        r'''Pretty-print available attributes to stderr.'''
        # global-var synchronization scheme only works with -j 1: adequate
        global RECORD_WAS_SHOWN
        if RECORD_WAS_SHOWN:
            return
        RECORD_WAS_SHOWN = True
        pp = pprint.PrettyPrinter(indent=4)
        # obs attributes
        if msg_obs is None:
            print('DRM observation attributes [use -a]:', file=sys.stderr)
        else:
            print(msg_obs, file=sys.stderr)
        # loop over keys => attributes appear as you'd name them with -a
        for k in sorted(obs.keys()):
            value_pp = pp.pformat(obs[k])
            print('  {}: {}'.format(k, value_pp), file=sys.stderr)
        if self.spc:
            print('', file=sys.stderr)
            if msg_spc is None:
                print('SPC attributes [use -S or -P]:', file=sys.stderr)
            else:
                print(msg_spc, file=sys.stderr)
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
            print('Note: SPC not loaded. Use --load_spc to force load.', file=sys.stderr)


    def extract_attrs(self, match, match_inv, attrs, show_attrs):
        r'''Extract attributes from each obs in the DRM.

        Returns a dictionary mapping attributes to lists, one list entry
        per DRM observation.
        '''
        # accumulate values in these lists, one per attribute
        vals = {name: [] for name in attrs.keys()}
        not_yet_shown = True
        for nobs, obs in enumerate(self.drm):
            # 1: no obs[match], or obs[match_inv] present => skip this observation
            if match and not self.extract_value_obs(obs, match, just_peeking=True):
                continue
            if match_inv and self.extract_value_obs(obs, match_inv, just_peeking=True):
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
            if VERBOSITY and Nplan > 1:
                print('----------')
                print(f'EX_ATTR: nobs = {nobs}, Nplan = {Nplan}')
                print(vals_1)
                #print(val_lens)
            # 3B: extend vals_1 along planets if needed
            # 3B.1 -- Nplan = 0 case
            #   skip this record if there was a planet attribute selected, but no planets there
            #   otherwise, we continue and the PlanetNum will tabulate as 0
            if Nplan == 0:
                if args.empty_skipped:
                    continue
                else:
                    # Nplan == 0 case, but not empty_skipped
                    # plan_num, if requested, will be 0
                    # ==> we will finish the loop, and output a record for this obs
                    #   fields will contain whatever extract_value() returned above
                    # fix up vals_1 for output
                    for name, attr in attrs.items():
                        # below: force all asked-for planet attributes to be [NaN]
                        # (at this point, they will be [])
                        # this singleton will be extend'ed onto vals below
                        if attr.flavor == 'planet':
                            vals_1[name] = [np.nan]
                        # further: correct a [] to [NaN]
                        elif len(vals_1[name]) == 0:
                            vals_1[name] = [np.nan]
                        # further: correct any [[]] to [NaN]
                        #   heuristic - planet attributes like SNR will be [[]]
                        elif len(vals_1[name]) == 1 and vals_1[name][0] == []:
                            vals_1[name] = [np.nan]
                    #print('Empty planets')
                    #print(vals_1)
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
        summary = self.extract_attrs(args.match, args.match_inv, args.all_attrs, args.show_attributes)
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
    sim = SimulationRun(f, load_spc=args.load_spc)
    #breakpoint()
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
            reductions = map_function(partial(outer_load_and_reduce, args=self.args, verb=self.args.verbose),
                                      self.sim_files)
        # hacky fix for no-drm case
        if len(self.sim_files) == 0:
            reductions = outer_load_and_reduce(None)
        # re-group the above reductions across sims
        # result is a dict containing reduced data, stored as self.summary
        self.regroup_and_accum(reductions, self.args.all_attrs.keys())

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
        r'''Dump reduced data to output.

        Args:
          args (namespace): program input arguments
          outfile (file): output filename
          '''
        saved_fields = args.all_attrs.keys()
        Nrows = set(len(self.summary[fname]) for fname in saved_fields)
        # Nrows will be a singleton set if all fields were in all DRMs
        # Nrows will be empty if no fields (e.g., -A)
        if len(Nrows) > 1:
            sys.stderr.write(f'Warning: DRMs have different attr counts: {sorted(Nrows)}\n')
        Nrow = max(Nrows, default=0)
        # if present, sort by seed to eliminate job-based ordering
        # sorted() is stable, so same-seed rows are not re-ordered
        inx = list(range(Nrow))
        if 'seed' in self.summary:
            # re-order the index list
            inx = sorted(inx,
                key=lambda i: self.summary['seed'][i])
        # make list-of-dicts to dump
        dumpable = []
        for i in inx:
            # dictionary mapping field -> value -- everything is a scalar here
            #print(f'{i} => {self.summary.seed}')
            #breakpoint()
            d = {f:self.summary[f][i] for f in saved_fields}
            dumpable.append(d)
        # place it in CSV or JSON
        if outfile:
            # default LW ~= 70, so small vectors would be line-wrapped in the CSV
            np.set_printoptions(linewidth=1000)
            if args.json:
                json.dump(dumpable, outfile, indent=2, cls=NumpyEncoder)
            else:
                csv_opts = {}
                if args.delimiter:
                    csv_opts['delimiter'] = args.delimiter
                w = csv.DictWriter(outfile, fieldnames=saved_fields, **csv_opts)
                if args.header:
                    w.writeheader()
                for d in dumpable:
                    for k,v in d.items():
                        # not used for ints or strings
                        if args.float_format and isinstance(v, (float, np.floating)):
                            d[k] = args.float_format % v
                    w.writerow(d)
        # list-of-dicts
        return dumpable


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
        sys.stderr.write(f'{args.progname}: Error loading json-format args from file "{args.file.name}"\n')
        sys.stderr.write(f'{args.progname}: Traceback follows.\n')
        raise
    # close file
    # remember only the name so "args" can be serialized by multiprocessing library
    args.file.close()
    args.file = args.file.name
    assert isinstance(d, collections.abc.Mapping), 'JSON "{}" must translate to a dictionary'.format(args.file)
    # 2: Propagate simple flags from d -> args
    # an args.match within given args overrides JSON
    if not args.match:
        args.match = d.get('__match__', args.match)
    # an args.match_inv within given args overrides JSON
    if not args.match_inv:
        args.match_inv = d.get('__match_inv__', args.match)
    # update args.header if present in d
    args.header = bool(d.get('__header__', args.header))
    # update args.load_spc if present in d
    args.load_spc = bool(d.get('__load_spc__', args.load_spc))
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
            name_colon_value = re.search(r'\s*(\w+[+]?)\s*:\s*(.*\S)\s*', argtext)
            if name_colon_value:
                # save the name, and the actual attr
                attr_name, attr_text = name_colon_value.groups()
            elif flavor in ('pseudo', 'attr', 'star', 'planet'):
                attr_name, attr_text = argtext, argtext.replace('+','')
            elif flavor == 'planets':
                if '+' in argtext:
                    raise ValueError(f'Vectorized attribute ({argtext}) illegal in "{flavor}"')
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
    # mapping, namely: flavor -> {name1:text1, name2:text2, ...}
    if args.file:
        # dictionary of the JSON file
        d = json_to_object(args)
        # make dict-of-dicts object
        obj = {}
        for flavor in AttrFlavors:
            obj[flavor] = d.get(f'__{flavor}__', {})
    else:
        # attributes-to-print are in argument-lists from the command line
        obj = argtexts_to_object(args)
    # Compose the list of all attributes we will extract (across all
    # flavors) by collecting attributes over flavors and names
    # Output field order preserved b/c dict is ordered (py 3.7+)
    all_attrs = dict()
    for flavor, attr_dict in obj.items():
        for name, attr_text in attr_dict.items():
            if name in all_attrs:
                sys.stderr.write(f'Warning: duplicate attribute {name}\n')
            p_expand = name.endswith('+')
            all_attrs[name] = Attribute(flavor, name.replace('+', ''), attr_text, p_expand)
    # the SPC must be loaded if any attributes requested it
    # (args.load_spc may already be True due to explicit CLI/JSON option)
    needs_spc = any(attr.flavor in ('star', 'planet', 'planets') for attr in all_attrs.values())
    if needs_spc:
        args.load_spc = True
    return all_attrs


########################################
###
###  main()
###
########################################

def main(args, infiles, outfile=None):
    # not going to allow for this case
    if len(infiles) == 0:
        sys.stderr.write(f'{args.progname}: No input files.\n')
        return {}
    if VERBOSITY:
        sys.stderr.write(f'{args.progname}: Loading {len(infiles)} files.\n')
    ensemble = EnsembleSummary(infiles, args)
    if VERBOSITY:
        sys.stderr.write(f'{args.progname}: Tabulating.\n')
    ensemble.load_and_reduce()
    d = ensemble.dump(args, outfile)
    return d


def parse_arglist(arglist):
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
    group.add_argument('-f', '--file', type=argparse.FileType('r'), help='JSON program file naming attributes')

    parser.add_argument('drm', metavar='DRM', nargs='*', default=[], help='drm file list')
    parser.add_argument('-m', '--match', type=str, default='', metavar='ATTR',
                            help='Attribute within obs to match, else, skip the obs')
    parser.add_argument('-M', '--match_inv', type=str, default='', metavar='ATTR',
                            help='Attribute within obs to NOT match (skip if present)')
    group = parser.add_argument_group('Output format')
    group.add_argument('-1', help='Supply CSV header line', action='store_true', dest='header',
                            default=False)
    group.add_argument('--delimiter', help='CSV field delimiter', metavar='STR', type=str)
    group.add_argument('--json', help='JSON output format', action='store_true', dest='json',
                            default=False)
    group = parser.add_argument_group('Seldom used')
    group.add_argument('--load_spc', help='Force load of SPC file', action='store_true', dest='load_spc',
                            default=False)
    group.add_argument('--format', help='Floating-point output format (%%-style)', dest='float_format', type=str, default='')
    group.add_argument('--empty', help='Output records despite having 0 planets', action='store_false', dest='empty_skipped',
                            default=True)
    group.add_argument('-A', help='Show list of Available Attributes on stderr', action='store_true', dest='show_attributes')
    group.add_argument('-v', help='Verbosity', action='count', dest='verbose', default=0)
    group.add_argument('-j', '--jobs', help='Number of parallel jobs (default = %(default)d)',
                      type=int, dest='jobs', default=int(N_CPU*0.65))
    # return the "args" namespace produced by parse_args()
    args = parser.parse_args(arglist)
    # for informative warnings, can be over-ridden
    args.progname = 'drm_tabulate'
    if args.float_format and '%' not in args.float_format:
        print(f'{args.progname}: Fatal. Need a % format specifier in --format.', file=sys.stderr)
        sys.exit(1)
    # Enforce that -P => --plan_num (only for CLI)
    if args.planet and 'plan_num' not in args.pseudo:
        args.pseudo.append('plan_num')
    # process all attr args into a unified list, args.all_attrs
    args.all_attrs = process_attr_program_inputs(args)
    return args


########################################
###
###  Program Entry
###
########################################

if __name__ == '__main__':
    args = parse_arglist(sys.argv[1:])
    
    # set umask in hopes that files will be group-writable
    os.umask(0o002)

    VERBOSITY = args.verbose

    # program name, for convenience
    args.progname = os.path.basename(sys.argv[0])
    
    # enforce that -A => -j 1 so stdout is not scrambled
    if args.show_attributes:
        args.jobs = 1

    if VERBOSITY > 0:
        sys.stderr.write('Attributes to be retrieved:\n')
        for name, attr in args.all_attrs.items():
            sys.stderr.write(f'  retrieving "{attr.flavor}" type attribute: {name} <- {attr.text}\n')

    main(args, args.drm, sys.stdout)
    sys.exit(0)
