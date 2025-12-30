#!/usr/bin/env python
r'''
reduce_drms.py: reduce a pile of DRMs to various summary CSV files

usage:
  `reduce_drms.py [ -O outfile ] [ -j N ] DRM [...]`

Args:
  DRM (file): a list of DRM pickles

Options:
  + `-O outfile` gives a template (containing exactly two occurrences of %s) for
     file outputs.  This is optional, but is best to supply.
  + `-j N` means to use N parallel workers to process the files.  By default,
     about 1/3 of the available cores will be used (25 on mustang*), which balances
     parallel makes (outside this routine) and parallel reductions (within it).
     If N = 0 or 1, no parallel workers are used, which is helpful for debugging.
  + `-D TYPE` gives a (string) TYPE that is inserted into the DEBUG variable in the 
     script, and which can be used to print a selected TYPE of debug output.
'''

# turmon jan 2018, oct 2018

## Note on design and code extensions
#
# This code is fundamentally a map/reduce over a list of DRMs:
#   for each DRM:
#     sim = load(SPC [i.e., star-planet configuration], DRM [i.e., observation-list])
#     summary = {}
#     summary.update(yield_analysis(sim))
#     [...on through many other kinds of analysis...]
#     summary.update(funnel_analysis(sim))
#     reductions.append(summary)
#   attrs = [list of keys in summary that we care about]
#   for attr in attrs:
#     accum[attr] = [r[attr] for r in reductions] # transpose
#   for attr in attrs:
#     summary[attr + '_mean'] = mean(accum[attr]) # <- unfortunate double-use of summary
#     [...on through _std, _q25, _q50, _q75, etc.]
#   dump_to_csv(summary)
#
# The initial "for" over DRMs is fast because it is parallelized and it is I/O heavy.
# The second "for" just transposes the array, and the last "for" uses standard numpy
# reductions (like np.mean, etc.) to boil down the list of numbers, one for each DRM,
# into a single summary.
#
# The single map/reduce principle embodied in this routine has worked well.  It allows 
# the DRM load (which is the most expensive part by far) to be done only once rather
# than spread over other executables.
#
# However, this file has grown larger than it should be.  It carries some technical
# debt and has needless repetition.
# Ways to retain the "single map/reduce" principle, but to disaggregate this
# file and improve DRY:
# (1) Establish a registry of X_analysis operations.  A driver would just run through
# the list, instead of manually running them one after the other.  This would simplify
# the summarize() method, and it would allow a uniform interface for each X_analysis()
# to be enforced.
# (2) In particular, the list of attributes ("attrs" above) is compiled semi-manually
# now.  This is a bad approach.  In some cases, we've used a (more) DRY approach of
# setting a special variable, named _X_keys for the X_analysis() method
# (e.g., _funnel_keys for funnel_analysis).  This variable holds the keys we care about.
# The registry above could build on this mechanism.
# (3) Namespace clashes within the summary[] dictionary are a potential problem.
# This could be addressed by not pooling the results of the X_analysis() routines
# into a single dictionary, it could be a dict-of-dicts, which might help
# modularity.
# (4) Similarly, the dump() method could be made more DRY, either by using the
# above idea of a "flavor of analysis that dumps a kind of file", or (more simply)
# a helper routine that dumps certain keys to a file of a given name.
#
# Note that one fundamental assumption is that the ensemble is a bag of interchangeable
# results.  We don't keep track of where a detection came from, and all the reductions
# (like mean, standard deviation, q50) are invariant to the ordering of the inputs.
#
# That said: To add another analysis mode to this goulash of analyses, do this:
# - Write a function, X_analysis, following delta_v_analysis as a model;
# - Put the keys it produces into a _X_keys key in the summary it returns;
# - Call the routine within summarize();
# - Add the "X" to the "key_family" in regroup_and_accum so that means, etc., are found;
# - Add the code segment to dump the summary from these keys to the dump() method.
# The relative ease of incrementally adding another analysis mode is of course the
# reason the above technical debt has accumulated.

import argparse
import sys
import shutil
import glob
import copy
import time
import json
import os
os.environ['OPENBLAS_NUM_THREADS'] = '4' # before "import numpy"
import gc
import csv
import re
import warnings
from functools import partial
from collections import defaultdict, Counter, namedtuple, OrderedDict
#import cPickle as pickle
import six.moves.cPickle as pickle
import six
#import multiprocessing.dummy
import multiprocessing as mproc
import numpy as np
from scipy import stats
from scipy.interpolate import interp1d
import astropy.units as u
import astropy.constants as const
from functools import reduce

########################################
##
## Globals
##
########################################

# global verbosity mode
VERBOSITY = 0
# global debug-output mode
DEBUG = ''

# number of CPUs
N_CPU = mproc.cpu_count()

# Temporal binning of mission elapsed time for detection-time plot
# (also for delta-v plot)
# default is 5-year mission, but see UPDATE_GLOBALS below
# np.arange(start_day, end_day, days_per_bin)
DETECTION_TIME_BINS = np.arange(0.0, 366*5.0, 30.5)

# Temporal binning of certain events (slews, chars) -- days
#   The histograms we compute for these events are normalized to sum to unity, 
#   so the bin-lists below have an extra bin to catch all the outliers, so that the
#   normalization is correct.  The extra bin is not output from this code.
DURATION_TIME_B0_BINS = np.concatenate((np.arange(0.0, 62.0, 2.0), [1e4]))
DURATION_TIME_B1_BINS = np.concatenate((np.arange(0.0, 31.0, 1.0), [1e4]))
DURATION_TIME_B2_BINS = np.concatenate((np.arange(0.0, 1.26, 1.0/24.0), [1e4]))
assert len(DURATION_TIME_B1_BINS) == len(DURATION_TIME_B2_BINS), 'Bins must be of same length (B1/B2)'
assert len(DURATION_TIME_B0_BINS) == len(DURATION_TIME_B1_BINS), 'Bins must be of same length (B0/B1)'

# Binning for Earth-chars
#   size is: "anticipated upper bound on number of Earth chars/mission"
EARTH_CHAR_COUNT_BINS = np.arange(0.0, 61.0, 1.0)

# Number of bins in slew/char/det histograms (bins are integers from 0..EVENT_COUNT_NBINS),
# size is controlled by "anticipated upper bound on number of detections/mission"
EVENT_COUNT_NBINS = 2000

# Temporal binning of promotion candidates
#  tailored for 5-year mission at <= 80% utilization
#  note, this is in units of instrument (detector) time
PROMOTION_TIME_BINS = np.arange(0.0, 366*5.0*0.80, 30.5)
# indexes (month numbers at current 30-day spacing) of T1 and T2
# These are in units of mission clock time
PROMOTION_PHIST_T1_INX = 24
PROMOTION_PHIST_T2_INX = 36
# promotion planet histograms go up to this many planets
PROMOTION_PHIST_NBINS = 40
PROMOTION_PHIST_BINS = np.arange(PROMOTION_PHIST_NBINS+1)

# filter incoming DRMs to characterization-only?  (disabled)
#MODE = 'char'
#MODE = 'det'
MODE = '*all*'

# bands we are collecting characterization info on
CHAR_BANDS = ['union', 'blue', 'red']


def UPDATE_GLOBALS(sim_info):
    r'''Update selected global variables upon initialization.'''
    global DETECTION_TIME_BINS
    missionLife = sim_info['missionLife']
    foo = copy.copy(DETECTION_TIME_BINS)
    DETECTION_TIME_BINS = np.arange(0.0, 366*missionLife, 30.5)


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

def np_force_string(x):
    r'''Convert np bytes_ array to a string array, if possible.

    Context: When Python3 reads a Python2 pickle with numpy strings, they will be
    bytes (i.e., encoded).  Before writing these strings (like star names) to CSV,
    we want them all to be python strings (but not unicode, which the csv-writer
    does not support).  We need to wrap in try/catch because for Python3 reading 
    Python3 pickles, these objects are already strings, and there is no decoding.
    '''
    try:
        return np.char.decode(x)
    except:
        return x

def char_within_band(char, band):
    r'''Does a particular characterization fall within the given notional band?'''
    if band == 'union':
        return True
    if (band == 'blue')  and (char['char_mode']['lam'] <= 600*u.nm):
        return True
    elif (band == 'red') and (char['char_mode']['lam'] >  600*u.nm):
        return True
    return False

def has_char_info(obs):
    r'''Utility function, does the observation have characterization info?

    See get_char_time for more on the need for this.'''
    return ('char_time' in obs) or ('char_info' in obs)

def get_char_time(obs):
    r'''Utility function, gets char time from a drm observation.

    The tricky issue here is that some characterizations (with starshade)
    have an embedded char_time key, and others (coronagraph + filters) have
    a list of char_info's with a series of char_times, which must be
    added up since they occur serially.'''
    if 'char_time' in obs:
        # starshade char - has a unitary char_time key
        return obs['char_time']
    else:
        # coronagraph-only char - sum integration windows
        return sum(cslice['char_time'] for cslice in obs['char_info'])

def get_char_status(obs):
    r'''Utility function, gets char status from a drm observation.

    Not yet used.'''
    if 'char_status' in obs:
        # starshade char - has a unitary char_status key
        # FIXME: is this the right place to find this key?  Is it not in char_info also?
        return obs['char_status']
    else:
        # coronagraph-only - get per-planet char status by combining channels
        # Each char status is a vector of -1, 0, +1 for each planet.
        # Planets are combined separately: 2 planets -> return a length-2 array.
        # There are 3x3 values to specify in combining 2 status values across
        # channels (wavelengths).  We combine according to these rules,
        # where s, s' are statuses {-1 (partial), 0 (fail), +1 (success)} --
        #   c(s,s) = s; c(s,s') = c(s',s); c(1, s) = 1; c(0, s) = s.
        # This turns out to be the same as the following:
        c = lambda s1, s2: np.sign(s1 + s2 + np.maximum(s1, s2))
        rv = 0.0 # start with the identity element
        for islice in obs['char_info']:
            rv = c(rv, islice['char_status'])
        return rv
    
# Fiddly helper function for JSON output
def array_encoder(obj):
    r"""Encodes numpy arrays, astropy Times, and astropy Quantities, into JSON.
    
    Called from json.dump for types that it does not already know how to represent,
    like astropy Quantity's, numpy arrays, etc.  The json.dump() method encodes types
    like integers, strings, and lists itself, so this code does not see these types.
    Likewise, this routine can and does return such objects, which is OK as long as 
    they unpack recursively into types for which encoding is known.
    
    """
    
    from astropy.time import Time
    from astropy.coordinates import SkyCoord
    if isinstance(obj, Time):
        # astropy Time -> time string
        return obj.fits # isot also makes sense here
    if isinstance(obj, u.quantity.Quantity):
        # note: it is possible to have a numpy ndarray wrapped in a Quantity.
        # NB: alternatively, can return (obj.value, obj.unit.name)
        return obj.value
    if isinstance(obj, SkyCoord):
        return dict(lon=obj.heliocentrictrueecliptic.lon.value,
                    lat=obj.heliocentrictrueecliptic.lat.value,
                    distance=obj.heliocentrictrueecliptic.distance.value)
    if isinstance(obj, (np.ndarray, np.number)):
        # ndarray -> list of numbers
        return obj.tolist()
    if isinstance(obj, (complex, np.complex)):
        # complex -> (real, imag) pair
        return [obj.real, obj.imag]
    if callable(obj):
        # this case occurs for interpolants like PSF and QE
        # We cannot simply "write" the function to JSON, so we make up a string
        # to keep from throwing an error.
        # The fix is simple: when generating the interpolant, add a _outspec attribute
        # to the function (or the lambda), containing (e.g.) the fits filename, or the
        # explicit number -- whatever string was used.  Then, here, check for that 
        # attribute and write it out instead of this dummy string.  (Attributes can
        # be transparently attached to python functions, even lambda's.)
        return 'interpolant_function'
    if isinstance(obj, set):
        return list(obj)
    if isinstance(obj, bytes):
        return obj.decode()
    # an EXOSIMS object
    if hasattr(obj, '_modtype'):
        return obj.__dict__
    # an object for which no encoding is defined yet
    #   as noted above, ordinary types (lists, ints, floats) do not take this path
    raise ValueError('Could not JSON-encode an object of type %s' % type(obj))


class RedirectStreams(object):
    r"""Set stdout and stderr to redirect to the named streams.

    Used for eliminating chatter to stdout upon module creation."""
    def __init__(self, stdout=None, stderr=None):
        self._stdout = stdout or sys.stdout
        self._stderr = stderr or sys.stderr

    def __enter__(self):
        self.old_stdout, self.old_stderr = sys.stdout, sys.stderr
        self.old_stdout.flush(); self.old_stderr.flush()
        sys.stdout, sys.stderr = self._stdout, self._stderr

    def __exit__(self, exc_type, exc_value, traceback):
        self._stdout.flush(); self._stderr.flush()
        sys.stdout = self.old_stdout
        sys.stderr = self.old_stderr

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
            
class YieldAccumulator(dict):
    r'''Custom accumulator with initial empty-list except for special cases.

    This accumulator allows us to accumulate various quantities into
    bins defined at runtime in yield_analysis().
    '''
    def __init__(self, *args, **kw):
        super(YieldAccumulator, self).__init__(*args, **kw)

    def __missing__(self, key):
        if type(key) is not str:
            raise KeyError(key) # disallow
        if key.startswith('exoE_'):
            # basic counter for exo-Earths
            new_value = 0.0
        elif key.startswith('set_'):
            # set accumulator for "is this a new observation"
            new_value = set()
        else:
            new_value = [] # default to list
        self[key] = new_value
        return new_value

    
class RpLBins(object):
    r'''Class to hold Rp - Luminosity bin properties.

    This is as much a container object for bin properties as it is a functional class.'''
    # Bin the detected planets into types
    # 1: planet-radius bin-edges  [units = Earth radii]
    # Old (early 2018, 3 bins):
    #   Rp_bins = np.array([0.5, 1.4, 4.0, 14.3])
    # New (May 2018, 5 bins x 3 bins, see Kopparapu et al, arxiv:1802.09602v1,
    # Table 1 and in particular Table 3 column 1, column 2 and Fig. 2):
    Rp_bins = np.array([0.5, 1.0, 1.75, 3.5, 6.0, 14.3])
    # 1b: bin lo/hi edges, same size as the resulting histograms
    # TODO: Luminosity should perhaps be in increasing order.  Plot how you want, but
    #       compute, store, and exchange in increasing order.
    Rp_lo = np.outer(Rp_bins[:-1], np.ones((3,1))).ravel()
    Rp_hi = np.outer(Rp_bins[1:],  np.ones((3,1))).ravel()
    # 2: stellar luminosity bins, in hot -> cold order
    #    NB: decreasing ordering.
    # Old (as above):
    # L_bins = np.array([
    #    [185, 1.5,  0.38, 0.0065],
    #    [185, 1.6,  0.42, 0.0065],
    #    [185, 1.55, 0.40, 0.0055]])
    # New (as above):
    L_bins = np.array([
        [182, 1.0,  0.28, 0.0035],
        [187, 1.12, 0.30, 0.0030],
        [188, 1.15, 0.32, 0.0030],
        [220, 1.65, 0.45, 0.0030],
        [220, 1.65, 0.40, 0.0025],
        ])
    # the below : selectors are correct for increasing ordering
    L_lo = L_bins[:,:-1]
    L_hi = L_bins[:,1:]
    # a bins: Unused
    #a_bins = [1./sqrt([185, 1.5, .38, .0065]),1./sqrt([185, 1.6, 0.42, .0065]),1./sqrt([185, 1.55, .4, .0055])]
    # total number of bins (e.g., 9 = 3*4-3, for a 3x3 histogram)
    RpL_bin_count = L_bins.size - (Rp_bins.size - 1)
    # radius/luminosity bin boundaries
    #   if there are 9 bins, there are 10 bin-edges, at 0.5, 1.5, ..., 9.5.
    #   this histogram drops the "0", or out-of-range, RpL region
    RpL_bin_edge_list = np.arange(0, RpL_bin_count+1) + 0.5

    # set up the bin-number map from (Rp_bin, L_bin) -> RpL_bin
    # this is a map from (int,int) -> int:
    #   yields 0 for input pairs outside the allowed range
    #        [namely, 1..len(Rp_bins) and 1..len(L_bins[i])]
    #   yields 1...9 otherwise, with 1, 2, 3 for the small planets (Rp bin number = 1).
    Rp_L_to_RpL_bin = defaultdict(int)
    # there are many ways to set this up: here is one.
    # (the below lines are enough for the old 3x3 setup, and they work for the new 5x3 setup too)
    Rp_L_to_RpL_bin[(1,1)] = 1 # smallest radius, highest luminosity => L bin number 1
    Rp_L_to_RpL_bin[(1,2)] = 2
    Rp_L_to_RpL_bin[(1,3)] = 3
    Rp_L_to_RpL_bin[(2,1)] = 4
    Rp_L_to_RpL_bin[(2,2)] = 5
    Rp_L_to_RpL_bin[(2,3)] = 6
    Rp_L_to_RpL_bin[(3,1)] = 7
    Rp_L_to_RpL_bin[(3,2)] = 8
    Rp_L_to_RpL_bin[(3,3)] = 9
    # New setup has 15 bins due to two new radius bins, so just add them here
    Rp_L_to_RpL_bin[(4,1)] = 10
    Rp_L_to_RpL_bin[(4,2)] = 11
    Rp_L_to_RpL_bin[(4,3)] = 12
    Rp_L_to_RpL_bin[(5,1)] = 13
    Rp_L_to_RpL_bin[(5,2)] = 14
    Rp_L_to_RpL_bin[(5,3)] = 15

    def quantize_orig(self, specs, plan_id, star_ind):
        r'''Compute the radius, luminosity, and combined bins for a given planet and star.

        This is Rhonda's original code. It is here to allow cross-checks but is not used now.
        Returns None if the planet/star lies outside the bin boundaries.'''

        # return value indicating error
        error_rval = (None, None, None)
        # extract planet values
        Rp_single = strip_units(specs['Rp'][plan_id])
        a_single = strip_units(specs['a'][plan_id])
        L_star = specs['L'][star_ind]
        L_plan = L_star/a_single**2
        # bin them
        tmpbinplace_Rp = np.digitize(Rp_single, self.Rp_bins)
        if tmpbinplace_Rp >= len(self.Rp_bins):
            return error_rval # out of range
        tmpbinplace_L = np.digitize(L_plan, self.L_bins[tmpbinplace_Rp-1])
        if tmpbinplace_L >= len(self.L_bins[tmpbinplace_Rp-1]):
            return error_rval # out of range
        Rp_bin = tmpbinplace_Rp.item()
        L_bin = tmpbinplace_L.item()
        RpL_bin = tmpbinplace_L.item() + 10*tmpbinplace_Rp.item()
        return Rp_bin, L_bin, RpL_bin

    def quantize(self, spc, plan_id, star_ind):
        r'''Compute the final radius/luminosity bin, an integer, for a given planet and star.

        Returns 0 if the planet/star lies outside the bin boundaries.  Returns 1..15 otherwise.
        Note: Not vectorized; scalar plan_id only.'''
        # extract planet and star properties
        Rp_plan = strip_units(spc['Rp'][plan_id])
        a_plan = strip_units(spc['a'][plan_id])
        L_star = spc['L'][star_ind]
        L_plan = L_star / (a_plan**2) # adjust star luminosity by distance^2 in AU
        # Bin by Rp.  "too low" maps to 0, and "too high" maps to len(Rp_bins).
        Rp_bin = np.digitize(Rp_plan, self.Rp_bins)
        # index into L_bins array: if Rp is out-of-range, index is irrelevant
        Rp_bin_index = Rp_bin-1 if (Rp_bin > 0) and (Rp_bin < len(self.Rp_bins)) else 0
        # bin by L
        L_bin = np.digitize(L_plan, self.L_bins[Rp_bin_index])
        # map the pair (Rp,L) -> RpL
        # Rp_bin and L_bin are np arrays, so need to cast to integers
        return self.Rp_L_to_RpL_bin[(int(Rp_bin), int(L_bin))]

    def is_hab_zone(self, spc, plan_id, star_ind):
        r'''Is the planet in the habitable zone?'''
        # Note: must work for bare integer plan_id, and for [] plan_id.
        try:
            if len(plan_id) == 0: return False
        except:
            pass # bare integer does not have len()
        # rescale SMA by luminosity
        L_star = spc['L'][star_ind]
        a_scaled = spc['a'][plan_id] / np.sqrt(L_star)
        return np.logical_and(
            a_scaled >=  .95*u.AU,
            a_scaled <= 1.67*u.AU)

    def is_earthlike(self, spc, plan_id, star_ind):
        r'''Is the planet earthlike?

        This version parallels the one in the EXOSIMS SurveySimulation prototype.'''
        try:
            if len(plan_id) == 0: return False
        except:
            pass # bare integer does not have len()
        scaleOrbits = True
        # extract planet and star properties
        Rp_plan = strip_units(spc['Rp'][plan_id])
        L_star = spc['L'][star_ind]
        if scaleOrbits:
            a_plan = strip_units(spc['a'][plan_id]) / np.sqrt(L_star)
        else:
            a_plan = strip_units(spc['a'][plan_id])
        # Definition: planet radius (in earth radii) and solar-equivalent luminosity must be
        # between the given bounds.
        Rp_plan_lo = 0.80/np.sqrt(a_plan)
        # We use the numpy versions so that plan_ind can be a numpy vector.
        return np.logical_and(
           np.logical_and(Rp_plan >= Rp_plan_lo, Rp_plan <= 1.4),
           np.logical_and(a_plan  >= 0.95,       a_plan  <= 1.67))

    def is_earthlike_old(self, spc, plan_id, star_ind):
        r'''Is the planet earthlike? -- Old version without orbit scaling.'''
        try:
            if len(plan_id) == 0: return False
        except:
            pass # bare integer does not have len()
        # extract planet and star properties
        Rp_plan = strip_units(spc['Rp'][plan_id])
        a_plan = strip_units(spc['a'][plan_id])
        L_star = spc['L'][star_ind]
        L_plan = L_star / (a_plan**2) # adjust star luminosity by distance^2 in AU
        # Definition: planet radius (in earth radii) and solar-equivalent luminosity must be
        # between the given bounds.
        # The magic numbers on L_plan are from:
        #    0.95 <= a/sqrt(L) <= 1.67 iff (1/1.67)^2 <= L/a^2 <= (1/0.95)^2
        # See also the condition in is_hab_zone, above.
        ## OLD:
        ## The lower Rp bound is not axis-parallel, but
        ## the best axis-parallel bound is 0.90, so that's what we use.
        ## Rp_plan_lo = 0.90
        # New: 0.8/sqrt(a)
        # NB: if scaleOrbits, the "a" here should instead be scaled by sqrt(L)!
        Rp_plan_lo = 0.80/np.sqrt(a_plan)
        # We use the numpy versions so that plan_ind can be a numpy vector.
        return np.logical_and(
            np.logical_and(Rp_plan >= Rp_plan_lo, Rp_plan <= 1.4),
            np.logical_and(L_plan  >= 0.3586,     L_plan  <= 1.1080))

class SimulationRun(object):
    r'''Load and summarize a simulation: one DRM and its corresponding SPC.'''
    def __init__(self, f, sim_info={}):
        # allow creating a dummy object so that its properties may be queried
        # (this is used in a couple of places)
        if f is None:
            self.drm = []
            self.spc = defaultdict(list)
            # set this one up (we query it later)
            self.spc['nPlans'] = 0
            self.name = 'dummy'
            self.Nstar = 0
            return
        # unpickling python2/numpy pickles within python3 requires this
        pickle_args = {} if sys.version_info.major < 3 else {'encoding': 'latin1'}
        # disabling gc during object construction speeds up by ~30% (12/2017, py 2.7.14)
        gc.disable()
        # drm = pickle.loads(open(f).read())
        drm = pickle.load(open(f, 'rb'), **pickle_args)
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
            # spc = pickle.loads(open(g).read())
            spc = pickle.load(open(g, 'rb'), **pickle_args)
        else:
            raise ValueError('Could not find a .spc file to match DRM <%s>' % f)
        # set up object state
        self.sim_info = sim_info # dict of information NOT from the DRM/SPC
        self.name = f
        self.Nstar = len(spc['Name'])
        self.spc = spc
        self.drm = drm
        self.summary = None # place-holder
    
    def summarize_revisits(self):
        r'''Compute revisit-count summary statistics (for one DRM).

        All of these results are binned by observation time.
        As a special case, returns the bin boundaries, which are presently just
        integer visit numbers in a range, if Nstar is given as 0
        Otherwise, returns the number of stars with the number of total visits
        falling into each visit-range bin.'''

        # maximum number of visits to histogram
        N_visit_max = 25
        # corresponding bin-boundaries -- includes upper boundary
        #   thus, visit_bins[-1] == N_visit_max
        visit_bins = np.arange(N_visit_max+1)
        # flag value to allow returning just the (lower) bin boundaries as a special case
        #   this is [0:N_visit_max-1]
        # Note: this is disabled -- we just put the bins in the return value in all cases
        if False and self.Nstar == 0:
            return {'h_visit_bins': visit_bins[:-1]}
        # for later, to determine if earthlike planet or not
        binner = RpLBins()
        # accumulate star indexes
        stars_visited_all   = []
        stars_visited_earth = []
        # DRM-FMT
        for obs in self.drm:
            # if it's a starshade DRM, and not a detection observation, then skip.
            if ('det_info' not in obs) and ('det_time' not in obs):
                continue 
            # record all stars-visited
            stars_visited_all.append(obs['star_ind'])
            # filter down to stars-with-any-earthlike
            if np.any(binner.is_earthlike(self.spc, np.array(obs['plan_inds']), obs['star_ind'])):
                stars_visited_earth.append(obs['star_ind'])
        # visit count, indexed by star (star index starts at 0)
        #  (replaced former np.int->np.int32 to silence warnings from np v1.19)
        visit_by_star_all   = np.bincount(np.array(stars_visited_all,   dtype=np.int32), minlength=self.Nstar)
        visit_by_star_earth = np.bincount(np.array(stars_visited_earth, dtype=np.int32), minlength=self.Nstar)
        # stars visited zero times, once, twice, etc., index from 0 to N_visit_max
        h_visit_all   = np.histogram(visit_by_star_all,   visit_bins)[0]
        h_visit_earth = np.histogram(visit_by_star_earth, visit_bins)[0]
        # adjust the final count to add in everything (strictly) above the top endpoint
        #  note, the top endpoint itself is already included in h_visit[-1]
        h_visit_all  [-1] += len(visit_by_star_all  [visit_by_star_all   > visit_bins[-1]])
        h_visit_earth[-1] += len(visit_by_star_earth[visit_by_star_earth > visit_bins[-1]])
        # return a dict of results
        rv = {
            'h_visit_all':   h_visit_all,
            'h_visit_earth': h_visit_earth,
            }
        rv['_summarize_revisits_keys'] = list(rv.keys())
        # statistics are not found over the bins: it's a trick to report binning
        # concept is to try to keep binning internal to this function if possible
        # TODO: report out these values differently!
        rv['h_visit_bins'] = visit_bins[:-1]
        return rv

    def per_star_yield(self):
        r'''Compute yield, tInt, etc., binned by target star (for one DRM).

        We currently account for detections (all), detections (Earth), 
        total integration time (detection, characterization),
        time-of-first observation (det/char), and a few other similar
        per-star summaries.
        Note: All of these summaries are itemized by star, not binned
        by radius and luminosity.
        '''

        # needed for is_earthlike()
        binner = RpLBins()
        # Accumulate counts into these variables - all indexed by star number
        # All are accumulators of various sorts - some work by addition and others by
        # maximum/minimum.
        # Implementation notes:
        # (1) Yield counts encompass the categories of:
        #   {det,char} X {plan,earth} X {cume,uniq}
        # (2) Yields are stored as vectors-of-vectors (VoV), and represent
        # *per-planet* yield counts for that star.  The VoV is a numpy 1xNstar
        # vector of numpy "Objects", each of which is a per-planet yield count.
        # (3) When a new set of detections is made at a star, we can do something like:
        #     yield[sind] += obs['det_status'], or
        #     yield[sind] = np.maximum(yield[sind], obs['det_status']),
        # and the per-planet detection counts will be updated.  By starting counts
        # at (scalar) zero, we don't have to special-case the first visit.
        # [A] detection = "det"
        n_star = self.Nstar
        dtime_ctr = np.zeros(n_star)
        tried_det_obs_time = np.nan + np.zeros(n_star) # time of first observation, or nan
        tried_det_ctr = np.zeros(n_star)
        det_comp = np.nan + np.zeros(n_star) # completeness, or nan
        # [A1] accumulators, across the DRM - for both cume and uniq
        yield_det_plan_ctr  = np.zeros(n_star, 'O') # VoV
        yield_det_earth_ctr = np.zeros(n_star, 'O') # VoV
        # [A2] Found after the DRM scan
        yield_det_plan_cume  = np.zeros(n_star)
        yield_det_plan_uniq  = np.zeros(n_star)
        yield_det_earth_cume = np.zeros(n_star)
        yield_det_earth_uniq = np.zeros(n_star)
        # [B] characterization = "char"
        ctime_ctr = np.zeros(n_star)
        tried_char_obs_time = np.nan + np.zeros(n_star) # time of first observation, or nan
        tried_char_ctr = np.zeros(n_star)
        char_comp = np.nan + np.zeros(n_star) # completeness, or nan
        # [B1] accumulators, across the DRM - for both cume and uniq
        yield_char_plan_ctr  = np.zeros(n_star, 'O') # VoV
        yield_char_earth_ctr = np.zeros(n_star, 'O') # VoV
        # [B2] Found after the DRM scan
        yield_char_plan_cume  = np.zeros(n_star)
        yield_char_plan_uniq  = np.zeros(n_star)
        yield_char_earth_cume = np.zeros(n_star)
        yield_char_earth_uniq = np.zeros(n_star)

        # DRM-FMT
        for obs in self.drm:
            sind = obs['star_ind']
            earths = binner.is_earthlike(self.spc, np.array(obs['plan_inds']), sind) # vector
            # catch detections in this clause
            if ('det_info' in obs) or ('det_time' in obs):
                # for coronagraph-only/Luvoir, detection info is kept in the 'det_info' list,
                # for starshade, detection info is in the drm entry itself (obs).
                # this abstracts the two cases by setting up a "pointer", obs_det.
                # but note, plan_inds and star_ind are always kept in obs itself.
                if 'det_info' in obs:
                    obs_det = obs['det_info'][0]
                else:
                    obs_det = obs
                # attempts-to-detect
                tried_det_ctr[sind] += 1
                # plug in this completeness (it may over-write an earlier one)
                if 'det_comp' in obs:
                    det_comp[sind] = obs['det_comp'] # no planets -> no comp
                # time of first detection attempt
                if tried_det_ctr[sind] == 1:
                    tried_det_obs_time[sind] = strip_units(obs['arrival_time'])
                # accumulate integration time, whether successful or not
                dtime_ctr[sind] += obs_det['det_time'].value
                # add 1 to the corresponding per-star planet count if detected
                yield_det_plan_ctr[sind] += (np.array(obs_det['det_status']) > 0)
                yield_det_earth_ctr[sind] += np.logical_and(earths, np.array(obs_det['det_status']) > 0)
            # catch characterizations in this clause
            if 'char_mode' in obs or 'char_info' in obs:
                # make "char_info" or a proxy of it ["char_info" is used in newer DRMs]
                #  char_info = [dict(char_time = X, char_status = Y) ...]
                if 'char_info' in obs:
                    char_info = obs['char_info']
                else:
                    char_info = [obs]
                # attempts-to-char
                tried_char_ctr[sind] += 1
                # plug in this completeness (it may over-write an earlier one)
                if 'char_comp' in obs:
                    char_comp[sind] = obs['char_comp'] # no planets -> no comp
                # time of first characterization attempt
                if tried_char_ctr[sind] == 1:
                    tried_char_obs_time[sind] = strip_units(obs['arrival_time'])
                # get just the first char in the list - don't double-count red + blue
                # TODO 6/2019: may need generalization to handle both coro-only + red/blue starshade cases?
                ## Formerly:
                ## ctime_ctr[sind] += strip_units(char_info[0]['char_time'])
                ctime_ctr[sind] += strip_units(get_char_time(obs))
                # find cumulative yield across multiple bands (for each planet)
                this_char_yield = np.zeros(1, 'bool') # i.e., false - loop will expand to vector
                for char in char_info:
                    char_status = np.array(char['char_status'])
                    # full-char (+1) or partial-char (-1) both count as a char
                    # multiple chars around one star are tracked separately
                    # accumulate here by "or"-ing all bands together
                    this_char_yield  = np.logical_or(this_char_yield, char_status != 0)
                # add the yield-vector (a 0/1 for each planet) to the per-star count vector
                yield_char_plan_ctr [sind] += this_char_yield
                yield_char_earth_ctr[sind] += np.logical_and(this_char_yield, earths)


        # summarize detections and characterizations, over whole DRM
        for sind in range(n_star):
            # sum *all* dets/chars of each planet around sind: [0, 1, 5, 2] -> 8
            yield_det_plan_cume [sind]  = np.sum(yield_det_plan_ctr[sind])
            yield_det_earth_cume[sind]  = np.sum(yield_det_earth_ctr[sind])
            yield_char_plan_cume [sind] = np.sum(yield_char_plan_ctr[sind])
            yield_char_earth_cume[sind] = np.sum(yield_char_earth_ctr[sind])
            # +1 for each planet around sind that has >0 dets/chars: [0, 1, 5, 2] -> 3
            yield_det_plan_uniq [sind]  = np.sum(yield_det_plan_ctr[sind]   > 0)
            yield_det_earth_uniq[sind]  = np.sum(yield_det_earth_ctr[sind]  > 0)
            yield_char_plan_uniq [sind] = np.sum(yield_char_plan_ctr[sind]  > 0)
            yield_char_earth_uniq[sind] = np.sum(yield_char_earth_ctr[sind] > 0)

        # count of #planets/star
        plan_per_star = 1.0*np.bincount(self.spc['plan2star'], minlength=n_star)
        earth_per_star = np.zeros(n_star)
        for sind in range(n_star):
            # all planets around star #sind
            plan_inds = np.where(self.spc['plan2star'] == sind)[0]
            # number of earthlike ones
            earth_per_star[sind] = np.sum(binner.is_earthlike(self.spc, plan_inds, sind))

        # "value" = (total yield) / (total time)
        # "frac" = (unique yield) / (# planets)
        # Suppress warnings on 0/0's here - NaN results are OK
        # In particular, NaNs in this DRM will be excluded from the overall average "value"
        # when we take the mean-across-DRMs later
        with np.errstate(divide='ignore', invalid='ignore'):
            det_plan_value   = yield_det_plan_cume   / dtime_ctr
            det_earth_value  = yield_det_earth_cume  / dtime_ctr
            char_plan_value  = yield_char_plan_cume  / ctime_ctr
            char_earth_value = yield_char_earth_cume / ctime_ctr
            det_plan_frac    = yield_det_plan_uniq   / plan_per_star
            det_earth_frac   = yield_det_earth_uniq  / earth_per_star
            char_plan_frac   = yield_char_plan_uniq  / plan_per_star
            char_earth_frac  = yield_char_earth_uniq / earth_per_star
            # turmon 2024/07: the average integration time, or NaN
            # (of course, != cumulative time on the target star)
            dtime_avg = dtime_ctr / tried_det_ctr
            ctime_avg = ctime_ctr / tried_char_ctr
        
        # return a dict of results for this DRM
        # Note: names from here forward do not change
        rv = {
            # det
            'h_star_det_visit':        tried_det_ctr,
            'h_star_det_tobs1':        tried_det_obs_time,
            'h_star_det_tInt':         dtime_ctr,
            'h_star_det_tIntAvg':      dtime_avg,
            'h_star_det_comp':         det_comp,
            'h_star_det_plan_cume':    yield_det_plan_cume,
            'h_star_det_plan_uniq':    yield_det_plan_uniq,
            'h_star_det_plan_value':   det_plan_value,
            'h_star_det_plan_frac':    det_plan_frac,
            'h_star_det_earth_cume':   yield_det_earth_cume,
            'h_star_det_earth_uniq':   yield_det_earth_uniq,
            'h_star_det_earth_value':  det_earth_value,
            'h_star_det_earth_frac':   det_earth_frac,
            # char
            'h_star_char_visit':       tried_char_ctr,
            'h_star_char_tobs1':       tried_char_obs_time,
            'h_star_char_tInt':        ctime_ctr,
            'h_star_char_tIntAvg':     ctime_avg,
            'h_star_char_comp':        char_comp,
            'h_star_char_plan_cume':   yield_char_plan_cume,
            'h_star_char_plan_uniq':   yield_char_plan_uniq,
            'h_star_char_plan_value':  char_plan_value,
            'h_star_char_plan_frac':   char_plan_frac,
            'h_star_char_earth_cume':  yield_char_earth_cume,
            'h_star_char_earth_uniq':  yield_char_earth_uniq,
            'h_star_char_earth_value': char_earth_value,
            'h_star_char_earth_frac':  char_earth_frac,
            # other
            'h_star_plan_per_star':    plan_per_star,
            'h_star_earth_per_star':   earth_per_star,
            }
        rv['_per_star_yield_keys'] = list(rv.keys())
        return rv

    def event_analysis(self):
        r'''Extract certain event information (slew, char, and det integration times) from DRM.

        Usage is binned by event duration.  This is distinct from resource_analysis(), which 
        bins according to mission elapsed time.
        The processing here is (presently) simple, and could have been folded in to another 
        loop over DRM events, but we preferred to factor this out for modularity.

        Algorithm:
        For each event class, keep a list of the duration of each event of that
        type across the DRM.  After scanning the full DRM, make a histogram of
        the durations.  This histogram will be averaged across the ensemble.
        '''

        # keep track of all event durations as a list
        event_det_duration = []
        event_char_duration = []
        event_slew_duration = []
        # DRM-FMT
        for obs in self.drm:
            # Process a detection
            if 'det_time' in obs or 'det_info' in obs:
                # for coronagraph-only/Luvoir, detection info is kept in the 'det_info' list,
                # for starshade, detection info is in the drm entry itself (obs).
                # this abstracts the two cases by setting up a "pointer", obs_det.
                # but note, plan_inds and star_ind are always kept in obs itself.
                if 'det_info' in obs:
                    obs_det = obs['det_info'][0]
                else:
                    obs_det = obs
                det_time = strip_units(obs_det['det_time'])
                event_det_duration.append(det_time)
            # Process a characterization
            if 'char_mode' in obs or 'char_info' in obs:
                if 'char_info' in obs:
                    char_info_1 = obs['char_info'][0]
                else:
                    char_info_1 = obs
                char_time = strip_units(get_char_time(obs)) # in days
                slew_time = strip_units(char_info_1.get('slew_time', 0.0)) # may not exist
                # skip "time = 0" chars: they are an artifact
                if char_time == 0: continue
                event_char_duration.append(char_time)
                event_slew_duration.append(slew_time)

        # bin the durations ("h_" is mnemonic for histogrammed)
        # We make densities out of these -- normalized so they integrate to unity.
        # The densities will be more useful than raw counts.  We have added, at the right endpoint,
        # a catch-all bin to contain very large durations so the true density will integrate to one.  
        # Then, for export, we chop off that last bin value (the trailing [:-1] below).
        # Empty lists (e.g., no slews) will throw a RuntimeWarning due to /0, which we suppress.
        # These lists will come through as all-NaN, and be eliminated from cross-ensemble
        # averages, which are computed with np.nanmean()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            h_event_char_b0_duration = np.histogram(event_char_duration, DURATION_TIME_B0_BINS, density=True)[0][:-1]
            h_event_slew_b0_duration = np.histogram(event_slew_duration, DURATION_TIME_B0_BINS, density=True)[0][:-1]
            h_event_det_b1_duration  = np.histogram(event_det_duration,  DURATION_TIME_B1_BINS, density=True)[0][:-1]
            h_event_char_b1_duration = np.histogram(event_char_duration, DURATION_TIME_B1_BINS, density=True)[0][:-1]
            h_event_slew_b1_duration = np.histogram(event_slew_duration, DURATION_TIME_B1_BINS, density=True)[0][:-1]
            h_event_det_b2_duration  = np.histogram(event_det_duration,  DURATION_TIME_B2_BINS, density=True)[0][:-1]
            h_event_char_b2_duration = np.histogram(event_char_duration, DURATION_TIME_B2_BINS, density=True)[0][:-1]
            h_event_slew_b2_duration = np.histogram(event_slew_duration, DURATION_TIME_B2_BINS, density=True)[0][:-1]

        # return value = dict with certain of the above variables
        namespace = locals()
        qois = [
            'h_event_char_b0_duration',
            'h_event_slew_b0_duration',
            'h_event_det_b1_duration',
            'h_event_char_b1_duration',
            'h_event_slew_b1_duration',
            'h_event_det_b2_duration',
            'h_event_char_b2_duration',
            'h_event_slew_b2_duration']
        rv = {qoi: namespace[qoi] for qoi in qois}
        # the list of keys that we're returning ... for use in later steps
        rv['_event_analysis_keys'] = qois
        return rv

    def count_events(self):
        r'''Extract event counts (#slew, #char, and #det) from DRM.

        The return value starts out as just a count (integer), and is then
        converted into a one-bin-on histogram, suitable for averaging
        across the ensemble.

        Algorithm:
        For each event class, keep a count of each event of that
        type across the DRM.  That count is then made into a one-bin-on
        histogram in the last lines.
        '''

        all_fields = ['det', 'char', 'slew', 'detp', 'det_rvplan', 'char_rvplan']
        # may not always be in the SPC dict
        promoted_stars = self.spc.get('promoted_stars', [])
        rv_plans = np.array(self.spc.get('known_earths', [])) # planet indexes
        # keep track of all event counts as a list
        events = Counter()
        # DRM-FMT
        for obs in self.drm:
            # Process a detection
            if 'det_time' in obs or 'det_info' in obs:
                events['det'] += 1
                if obs['star_ind'] in promoted_stars:
                    events['detp'] += 1
                # we don't care if more than 1 planet, or success/fail
                p_detected = np.array(obs['plan_inds'])
                if len(np.intersect1d(rv_plans, p_detected)) > 0:
                    events['det_rvplan'] += 1
            # Process a characterization
            if 'char_mode' in obs or 'char_info' in obs:
                p_detected = np.array(obs['plan_inds'])
                if 'char_info' in obs:
                    char_info_1 = obs['char_info'][0]
                else:
                    char_info_1 = obs
                char_time = strip_units(get_char_time(obs))
                slew_time = strip_units(char_info_1.get('slew_time', 0.0)) # may not exist
                # skip "time = 0" chars: they are an artifact
                if char_time > 0:
                    events['char'] += 1
                    if len(np.intersect1d(rv_plans, p_detected)) > 0:
                        events['char_rvplan'] += 1
                # if slew_time was not there, it will have been set to 0 above
                if slew_time > 0:
                    events['slew'] += 1
        # make a dict of 0/1 histograms, which is convenient to average over DRMs later
        final_fields = []
        rv = {}
        for event_name in all_fields:
            field = 'h_event_count_%s' % event_name
            final_fields.append(field)
            rv[field] = np.zeros(EVENT_COUNT_NBINS)
            rv[field][min(events[event_name], EVENT_COUNT_NBINS-1)] = 1.0
        # the list of keys that we're returning ... for use in later steps
        rv['_event_count_keys'] = final_fields
        # return value = above dict
        return rv


    def resource_analysis(self):
        r'''Extract time-binned resource usage (fuel and integration time) from DRM.

        Usage is binned by mission elapsed time.  This routine is written so that 
        long-duration events, with resource usage linear over the event, have the usage
        spread out evenly over the duration of the event.  This includes fuel usage and 
        integration time. 
        In progress:
        We also keep track of detections, but we do not spread these out over time.'''

        # handle the coronagraph-only case, where fuel use is not in the DRM
        # FIXME: Since this routine was generalized from "fuel" to "resources",
        # we can't just quick-exit on no-fuel-present like this.
        # For now, the fix is to eliminate this IF and suffer the consequences,
        # (if any!) and then fix the main code, as needed, for current DRMs.
        if False:
            zero = np.zeros(len(DETECTION_TIME_BINS))
            return {
                'h_time_det_cume':  zero,
                'h_time_det_incr':  zero,
                'h_time_char_cume': zero,
                'h_time_char_incr': zero,
                'h_time_slew_cume': zero,
                'h_time_slew_incr': zero,
                'h_time_fuel_slew': zero,
                'h_time_fuel_keep': zero,
                'h_time_fuel_all':  zero,
                }

        mismatch = 0
        n_slew = 0
        # initialize station-keeping:
        #   (sk_time, sk_fuel) earliest point = (-1,0)
        #   (sk_time, char_time_cume) earliest point = (-1,0)
        sk_time_tiepoint = [-1.0]
        sk_fuel_tiepoint = [ 0.0]
        char_time_cume_tiepoint = [0.0]
        # initialize slew:
        #   (slew_time, slew_fuel) earliest point = (-1,0)
        #   (slew_time, slew_time_cume) earliest point = (-1,0)
        slew_time_tiepoint = [-1.0]
        slew_fuel_tiepoint = [ 0.0]
        slew_time_cume_tiepoint = [0.0]
        # initialize detections:
        #   (det_time, det_time_cume) earliest point = (-1,0)
        det_time_tiepoint = [-1.0]
        det_time_cume_tiepoint = [0.0]

        # Algorithm:
        # For fuel usage for both slew and station-keeping, make a list of tie points of
        # cumulative fuel usage.  Then interpolate linearly between these tie points.
        # For SK, the tie points are (arrival_time, arrival_time + char_time).
        # For slew, tie points are  (arrival_time - slew_time, arrival_time)
        # Fuel is expended linearly in the tie point interval, and the value for the
        # left end of the tie point is taken from the prior value.
        # Obtain monthly usage by making a linear interpolator governed by these time
        # points.  This allows a consistent "time axis" for averaging across DRMs.

        # a small time (in days), to prevent instantaneous jumps in fuel use
        tiny_time = 1e-6
        # DRM-FMT
        for obs in self.drm:
            # DRM-FMT
            arrival_time = strip_units(obs['arrival_time'])
            if 'det_time' in obs or 'det_info' in obs:
                # for coronagraph-only/Luvoir, detection info is kept in the 'det_info' list,
                # for starshade, detection info is in the drm entry itself (obs).
                # this abstracts the two cases by setting up a "pointer", obs_det.
                # but note, plan_inds and star_ind are always kept in obs itself.
                if 'det_info' in obs:
                    obs_det = obs['det_info'][0]
                else:
                    obs_det = obs
                det_time_1 = strip_units(obs_det['det_time'])
                prior_time = det_time_tiepoint[-1]
                assert arrival_time >= prior_time, 'detection time sequence mismatch'
                det_time_0 = det_time_cume_tiepoint[-1]
                det_time_tiepoint.extend([arrival_time, arrival_time + max(tiny_time,det_time_1)])
                det_time_cume_tiepoint.extend([det_time_0, det_time_0 + det_time_1])

            # Process a characterization
            if 'char_mode' in obs or 'char_info' in obs:
                # for multi-mode char (Luvoir), obs may be a list - take the first
                if 'char_info' in obs:
                    char_info_1 = obs['char_info'][0]
                else:
                    char_info_1 = obs
                n_slew = n_slew + 1
                slew_time = strip_units(char_info_1.get('slew_time', 0.0))
                char_time = strip_units(get_char_time(obs)) # in days
                # fuel use for slew
                #   ramp from arrival_time-slew_time to arrival_time
                #   max(...) ensures no instantaneous jump, because slew_time = 0 at start
                prior_time = slew_time_tiepoint[-1]
                slew_time_cume_0 = slew_time_cume_tiepoint[-1]
                slew_fuel_0 = slew_fuel_tiepoint[-1]
                slew_fuel = strip_units(char_info_1.get('slew_mass_used', 0.0))
                slew_time_tiepoint.extend([max(arrival_time - slew_time, prior_time+tiny_time), arrival_time])
                slew_fuel_tiepoint.extend([slew_fuel_0, slew_fuel_0 + slew_fuel])
                slew_time_cume_tiepoint.extend([slew_time_cume_0, slew_time_cume_0 + slew_time])
                # assert that this pair of tie points lies after the prior one
                # FIXME: should not happen (?)
                if (arrival_time - slew_time) < prior_time:
                    mismatch += 1
                    #import pdb; pdb.set_trace()
                #assert (arrival_time - slew_time) >= prior_time, 'slew time sequence mismatch'
                # fuel use for station-keeping
                #   ramp over course of observation
                #   max(...) ensures no instantaneous jump
                prior_time = sk_time_tiepoint[-1]
                #assert arrival_time >= prior_time, 'station-keeping time sequence mismatch'
                if False and arrival_time < prior_time:
                    print('Stopping due to stationkeeping anomaly')
                    import pdb; pdb.set_trace()
                char_time_0 = char_time_cume_tiepoint[-1]
                sk_fuel_0 = sk_fuel_tiepoint[-1]
                # if the char_mass is not recorded, just use 0
                # (we keep track of char_time_cume + char_mass using the same tiepoints)
                sk_fuel = strip_units(char_info_1.get('char_mass_used', 0.0))
                sk_time_tiepoint.extend([arrival_time, arrival_time + max(tiny_time,char_time)])
                sk_fuel_tiepoint.extend([sk_fuel_0, sk_fuel_0 + sk_fuel])
                char_time_cume_tiepoint.extend([char_time_0, char_time_0 + char_time])

        if False and mismatch > 0:
            print('Observed %d unusual slews' % mismatch)
        # account for empty DRMs by adding a final tiepoint
        if len(slew_time_tiepoint) == 1:
            slew_time_tiepoint.append(DETECTION_TIME_BINS[-1])
            slew_fuel_tiepoint.append(0.0)
            slew_time_cume_tiepoint.append(0.0)
        if len(sk_time_tiepoint) == 1:
            sk_time_tiepoint.append(DETECTION_TIME_BINS[-1])
            sk_fuel_tiepoint.append(0.0)
            char_time_cume_tiepoint.append(0.0)
        if len(det_time_tiepoint) == 1:
            det_time_tiepoint.append(DETECTION_TIME_BINS[-1])
            det_time_cume_tiepoint.append(0.0)
        # extend to end of mission timeline - py 2.7.14 interp1d does not support fill values consistently
        final_time = DETECTION_TIME_BINS[-1]
        det_time_tiepoint.append(final_time) # i1
        slew_time_tiepoint.append(final_time) # i2, i4
        sk_time_tiepoint.append(final_time) # i3, i5
        det_time_cume_tiepoint.append(det_time_cume_tiepoint[-1]) # i1
        slew_time_cume_tiepoint.append(slew_time_cume_tiepoint[-1]) # i2
        char_time_cume_tiepoint.append(char_time_cume_tiepoint[-1]) # i3
        slew_fuel_tiepoint.append(slew_fuel_tiepoint[-1]) # i4
        sk_fuel_tiepoint.append(sk_fuel_tiepoint[-1]) # i5
        # i1: linear interpolator for time-spent-detecting
        det_time_interp_func = interp1d(det_time_tiepoint, det_time_cume_tiepoint)
        h_time_det = det_time_interp_func(DETECTION_TIME_BINS)
        # i2: linear interpolator for time-spent-slewing
        slew_time_interp_func = interp1d(slew_time_tiepoint, slew_time_cume_tiepoint)
        h_time_slew = slew_time_interp_func(DETECTION_TIME_BINS)
        # i3: linear interpolator for time-spent-characterizing
        char_time_interp_func = interp1d(sk_time_tiepoint, char_time_cume_tiepoint)
        h_time_char = char_time_interp_func(DETECTION_TIME_BINS)
        # i4: linear interpolator for slew fuel
        slew_fuel_interp_func = interp1d(slew_time_tiepoint, slew_fuel_tiepoint)
        h_time_fuel_slew = slew_fuel_interp_func(DETECTION_TIME_BINS)
        # i5: linear interpolator for station-keeping fuel
        sk_fuel_interp_func = interp1d(sk_time_tiepoint, sk_fuel_tiepoint)
        h_time_fuel_keep = sk_fuel_interp_func(DETECTION_TIME_BINS)
        # combined fuel use
        h_time_fuel_all = h_time_fuel_slew + h_time_fuel_keep
        # return value
        rv = {
            'h_time_det_cume':  h_time_det,
            'h_time_det_incr':  np.ediff1d(h_time_det, to_begin=0.0),
            'h_time_char_cume': h_time_char,
            'h_time_char_incr': np.ediff1d(h_time_char, to_begin=0.0),
            'h_time_slew_cume': h_time_slew,
            'h_time_slew_incr': np.ediff1d(h_time_slew, to_begin=0.0),
            'h_time_fuel_slew': h_time_fuel_slew,
            'h_time_fuel_keep': h_time_fuel_keep,
            'h_time_fuel_all':  h_time_fuel_all,
            }
        rv['_resource_analysis_keys'] = list(rv.keys())
        return rv

    def delta_v_analysis(self):
        r'''Extract time-binned delta-V information.

        This analysis is handled differently than the earlier resource_analysis() algorithm.
        That setup used a linear ramp of fuel expenditure across the observation window,
        which proved difficult to manage because sometimes windows could overlap or be
        very narrow.  So here, we parcel delta-v into a discrete number of chunks that 
        associated with an observation window, and later bin those chunks into a histogram.
        '''

        def parcel_delta_v(t0, dt, dv):
            r'''Split dv into a given number of parcels, spread over "dt" time units.'''
            dv_m_s = strip_units(dv)
            # let the #parcels scale with dt -- figuring 2-day bins is the most we need
            N_parcel = max(3, int(np.ceil(dt / 2.0)))
            # parcel out time and dv
            t_parcel = np.linspace(t0, t0+max(0.0, dt), N_parcel)
            dv_parcel = dv_m_s/float(N_parcel) + np.zeros(t_parcel.shape)
            return t_parcel, dv_parcel

        # compose lists of times and corresponding delta-v
        slew_times, slew_dvs = [], []
        det_times,  det_dvs  = [], []
        char_times, char_dvs = [], []
        for obs in self.drm:
            # DRM-FMT
            arrival_time = strip_units(obs['arrival_time'])
            if 'slew_dV' in obs and 'slew_time' in obs:
                slew_time = strip_units(obs['slew_time'])
                # TODO: should be parceled over (arrival_time-slew_time, arrival_time)
                t1, dv1 = parcel_delta_v(arrival_time, slew_time, obs['slew_dV'])
                slew_times.extend(t1)
                slew_dvs.extend(dv1)
            if 'det_dV' in obs and 'det_time' in obs:
                det_time = strip_units(obs['det_time'])
                t1, dv1 = parcel_delta_v(arrival_time, det_time, obs['det_dV'])
                det_times.extend(t1)
                det_dvs.extend(dv1)
            if 'char_dV' in obs:
                char_time = strip_units(get_char_time(obs)) # char_time can be in a list under obs
                t1, dv1 = parcel_delta_v(arrival_time, char_time, obs['char_dV'])
                char_times.extend(t1)
                char_dvs.extend(dv1)

        # these histograms are all incremental (per-month)
        # ...binned by time, weighted by delta-v
        h_time_delta_v_slew = np.histogram(slew_times, bins=DETECTION_TIME_BINS, weights=slew_dvs)[0]
        h_time_delta_v_det  = np.histogram(det_times,  bins=DETECTION_TIME_BINS, weights=det_dvs) [0]
        h_time_delta_v_char = np.histogram(char_times, bins=DETECTION_TIME_BINS, weights=char_dvs)[0]
        # combined delta-v for observations taken by the telescope bus
        h_time_delta_v_obs = h_time_delta_v_det + h_time_delta_v_char

        # list of all histograms we want to propagate outward
        all_keys = ['h_time_delta_v_%s' % (key, ) for key in ('slew', 'det', 'char', 'obs')]
        # compose a dictionary holding result
        #   ...maintain key order so we can extract its keys later
        rv = OrderedDict()
        sources = locals()
        for key in all_keys:
            # histogram of delta-v, binned by time
            h1 = sources[key]
            # save this "incremental" histogram and its cumulative version
            rv['%s_incr' % key] = h1
            rv['%s_cume' % key] = np.cumsum(h1)
        # return value
        rv_keys = list(rv.keys()) # we're about to add one key, and we don't want it in this list
        rv['_delta_v_keys'] = rv_keys
        return rv

    def funnel_analysis(self):
        r'''Extract promotion and characterization counts from DRM.

        For each sim, these are scalar quantities.'''
        ## 0: set up return value
        rv = dict()
        # initialize counters -- force them to exist
        for list_name in ('promo', 'deep'):
            for target in ('star', 'allplan', 'hzone', 'earth'):
                rv['funnel_%s_%s' % (list_name, target)] = 0.0
                for outcome in ('tries', 'chars'):
                    for count in ('cume', 'uniq'):
                        rv['funnel_%s_%s_%s_%s' % (list_name, outcome, target, count)] = 0.0
        # record all the above keys for later use -- needs to be an explicit list
        rv['_funnel_keys'] = list(rv.keys())
        # add in the list of earth char attempts
        rv['earth_char_list'] = []
        rv['_earth_char_keys'] = ['earth_char_list']
        # Do not error-out if sim load did not work, resulting in sim_info being not present or None
        if not getattr(self, 'sim_info', False):
            return rv
        binner = RpLBins() # for is_earthlike()
        # earthlike = binner.is_earthlike(self.spc, np.arange(self.spc['nPlans']), self.spc['plan2star'])
        # hzone = binner.is_hab_zone(self.spc, np.arange(self.spc['nPlans']), self.spc['plan2star'])

        # Apparent magnitude of stars (from info in the .spc)
        #  4.83 is the absolute visual magnitude of the Sun (mag for L = 1, d = 10pc)
        #  The other two terms correct for luminosity and distance
        #  L is in solar luminosity units (L increases -> Vmag decreases)
        #  distance: 2.5*log10(dist^2/10pc) = 5 * [log10(dist) - 1]  (dist incr. -> Vmag incr.)
        #  (we do not have a bolometric correction for each star, so by using 4.83, we are in effect
        #  using the Sun's bolometric correction for all stars)
        Vmag = 4.83 - 2.5 * np.log10(self.spc['L']) + 5.0 * (np.log10(self.spc['dist'].to('pc').value) - 1.0)

        # Up here so it's available throughout
        # 'promoted_stars' may not always be in the SPC dict
        promoted_stars = self.spc.get('promoted_stars', [])
        # deep-dive stars
        top_HIPs = self.sim_info['top_HIPs'] # list of hipparcos names of deep-dive stars
        top_sInds = np.where(np.isin(self.spc['Name'], top_HIPs))[0]

        ## 1: Find char yield at each star, for each planet, across the DRM
        # -- Some yields are stored as vectors-of-vectors (VoV), and represent
        # *per-planet* yield counts for that star.  The VoV is a numpy 1xNstar
        # vector of numpy "Objects", each of which is a per-planet yield count.
        # -- When a new set of chars (1xNplan) is made at a star, we:
        #     yield[sind] += this_char_yield
        # and the per-planet char counts will be updated.  By starting counts
        # at (scalar) zero, we don't have to special-case the first visit.
        n_star = self.Nstar
        tried_char_ctr = np.zeros(n_star)
        yield_char_ctr = np.zeros(n_star)
        yield_char_plan_ctr  = np.zeros(n_star, 'O') # VoV
        yield_char_hzone_ctr = np.zeros(n_star, 'O') # VoV
        yield_char_earth_ctr = np.zeros(n_star, 'O') # VoV
        for obs in self.drm:
            # skip detections
            #if ('det_info' in obs) or ('det_time' in obs):
            if not has_char_info(obs):
                continue
            # skip chars that were cancelled
            if get_char_time(obs) == 0.0:
                continue
            # (only chars below here)
            sind = obs['star_ind']
            plan_inds = np.array(obs['plan_inds']) # vector
            earths = binner.is_earthlike(self.spc, plan_inds, sind) # vector
            hzones = binner.is_hab_zone( self.spc, plan_inds, sind) # vector

            if ('yield' in DEBUG) and (sind in top_sInds):
                print('Observed a DD target %3d -- #EEC = %d -- chartime = %g ' % (
                    sind, np.sum(earths), get_char_time(obs).value))

            # make "char_info" or a proxy of it ["char_info" is used in newer DRMs]
            #  char_info = [dict(char_time = X, char_status = Y) ...]
            if 'char_info' in obs:
                char_info = obs['char_info']
            else:
                char_info = [obs]
            # attempts-to-char
            tried_char_ctr[sind] += 1
            # find cumulative yield across multiple bands (for each planet)
            this_char_yield = np.zeros(1, 'bool') # i.e., false - loop will expand to vector
            for char in char_info:
                char_status = np.array(char['char_status'])
                # full-char (+1) or partial-char (-1) both count as a char
                # multiple chars around one star are tracked separately
                # accumulate here by "or"-ing all bands together
                this_char_yield  = np.logical_or(this_char_yield, char_status != 0)
            # adds 1 if any planet had a successful char
            yield_char_ctr[sind] += np.any(this_char_yield)
            # add the yield-vector (a 0/1 for each planet) to the per-star count vector
            yield_char_plan_ctr [sind] += this_char_yield
            yield_char_hzone_ctr[sind] += np.logical_and(this_char_yield, hzones)
            yield_char_earth_ctr[sind] += np.logical_and(this_char_yield, earths)
            # report on earth characterization attempts
            # short-circuit this if no earths
            if np.any(earths):
                sind_promo = sind in promoted_stars
                sind_deep  = sind in top_sInds
                # NB: obs[char_params][WA] != obs[char_WA] (5/2019: believe fixed in Exosims?)
                parms = obs['char_params'] if 'char_params' in obs else char_info[0]['char_params']
                earth_inxs = np.where(earths)[0]
                # one output record per earth
                for earth_inx in earth_inxs:
                    # FIXME: is there a better way than path manipulation to get the ensemble number?
                    try:
                        ensemble_num = int(os.path.splitext(os.path.basename(self.name))[0])
                    except:
                        ensemble_num = -1
                    if 'char_SNR' in obs:
                        char_SNR_1 = obs['char_SNR'][earth_inx]
                    else:
                        char_SNR_1 = char_info[0]['char_SNR'][earth_inx]
                    # spc['Name'][] is a bytes sequence: decode into string for output
                    char_dict = OrderedDict([
                        ('ensemble', ensemble_num),
                        ('obsnum', obs['Obs#'] if 'Obs#' in obs else obs['ObsNum']),
                        ('name', np_force_string(self.spc['Name'][sind])), # .decode('utf-8')
                        ('sind', sind),
                        ('pind', earth_inx),
                        ('n_earth',    int(np.sum(earths))),
                        ('n_success',  int(np.sum(this_char_yield[earth_inxs]))),
                        ('is_success', int(this_char_yield[earth_inx])),
                        ('is_deep', int(sind_deep)),
                        ('is_promo', int(sind_promo)),
                        ('WA',   parms['WA']    [earth_inx].to('mas').value),
                        ('dMag', parms['dMag']  [earth_inx]),
                        ('phi',  parms['phi']   [earth_inx]),
                        ('char_SNR', char_SNR_1),
                        ('MV', Vmag[sind]),
                        ])
                    rv['earth_char_list'].append(char_dict)
                    # FIXME: temporary, for investigating char fails (2/2019)
                    if 'yield' in DEBUG:
                        if (sind in top_sInds):
                            print('%s,%s,%d,%d,%d,%d,%f,%f,%f,%f' % (self.spc['Name'][sind],
                                                    os.path.splitext(os.path.basename(self.name))[0],
                                                    obs['Obs#'] if 'Obs#' in obs else obs['ObsNum'],
                                                    1 if this_char_yield[earth_inx] else 0,
                                                    sind, earth_inx,
                                                    parms['WA']    [earth_inx].to('mas').value,
                                                    parms['dMag']  [earth_inx],
                                                    parms['phi']   [earth_inx],
                                                    char_SNR_1,
                                                          ))
        
        ## 2: Use above yields to compile metrics for promoted stars
        for star_list, list_name in ((promoted_stars, 'promo'), (top_sInds, 'deep')):
            for sind in star_list:
                # put the list_name into the template where indicated
                def n(template):
                    return template % list_name
                # planets for star #sind
                plan_inds = np.where(self.spc['plan2star'] == sind)[0]
                hzones = binner.is_hab_zone (self.spc, plan_inds, sind)
                earths = binner.is_earthlike(self.spc, plan_inds, sind)

                ## FIXME: temporary, for testing deep-dive chars (2/2019)
                if False and list_name == 'deep':
                    print(",".join((os.path.basename(self.name),
                                        self.spc['Name'][sind], '%d' % np.sum(earths), '%d' % len(plan_inds))))

                ## 2A: promotion counts
                rv[n('funnel_%s_star')]    += 1.0
                rv[n('funnel_%s_allplan')] += 1.0 * len(plan_inds)
                rv[n('funnel_%s_hzone')]   += 1.0 * np.sum(hzones)
                rv[n('funnel_%s_earth')]   += 1.0 * np.sum(earths)

                ## 2B: char attempts
                # char attempt/star
                rv[n('funnel_%s_tries_star_cume')] += tried_char_ctr[sind]
                rv[n('funnel_%s_tries_star_uniq')] += np.minimum(1, tried_char_ctr[sind])
                # char attempt/allplan
                rv[n('funnel_%s_tries_allplan_cume')] += tried_char_ctr[sind] * len(plan_inds)
                rv[n('funnel_%s_tries_allplan_uniq')] += np.minimum(1, tried_char_ctr[sind]) * len(plan_inds)
                # char attempt/hzone
                rv[n('funnel_%s_tries_hzone_cume')] += tried_char_ctr[sind] * np.sum(hzones)
                rv[n('funnel_%s_tries_hzone_uniq')] += np.minimum(1, tried_char_ctr[sind]) * np.sum(hzones)
                # char attempt/earth
                #   (if 2 earths X 3 tries, uniq = 2; unique means at least one success across tries)
                rv[n('funnel_%s_tries_earth_cume')] += tried_char_ctr[sind] * np.sum(earths)
                rv[n('funnel_%s_tries_earth_uniq')] += np.minimum(1, tried_char_ctr[sind]) * np.sum(earths)

                ## 2B: char success
                # char success/star
                rv[n('funnel_%s_chars_star_cume')]  += yield_char_ctr[sind]
                rv[n('funnel_%s_chars_star_uniq')]  += np.minimum(1, yield_char_ctr[sind])
                # char success/allplan
                rv[n('funnel_%s_chars_allplan_cume')] += np.sum(yield_char_plan_ctr[sind])
                rv[n('funnel_%s_chars_allplan_uniq')] += np.sum(yield_char_plan_ctr[sind] > 0)
                # char success/hzone
                rv[n('funnel_%s_chars_hzone_cume')] += np.sum(yield_char_hzone_ctr[sind])
                rv[n('funnel_%s_chars_hzone_uniq')] += np.sum(yield_char_hzone_ctr[sind] > 0)
                # char success/earth
                rv[n('funnel_%s_chars_earth_cume')] += np.sum(yield_char_earth_ctr[sind])
                rv[n('funnel_%s_chars_earth_uniq')] += np.sum(yield_char_earth_ctr[sind] > 0)

        ## 3: return the dict of results for this DRM
        # Note: valid keys are in '_funnel_keys'
        return rv


    def det_funnel_analysis(self):
        r'''Extract detection observation status from DRM.

        For each sim, these are scalar quantities.'''
        ## 0: set up return value
        rv = dict()
        # target types - star is a bit special but can be handled similarly
        Targets = ('star', 'allplan', 'hzone', 'earth')
        # initialize counters -- force them to exist
        # FIXME: rework to use f-strings
        for sname in ('allstar', 'promo'):
            for target in Targets:
                # tries, success
                rv[f'detfunnel_{sname}_cand_{target}'] = 0.0
                for outcome in ('tries', 'fails', 'success'):
                    for count in ('cume', 'uniq'):
                        rv[f'detfunnel_{sname}_{outcome}_{target}_{count}'] = 0.0
                count = 'cume'
                for outcome in ('f_snr', 'f_iwa', 'f_owa'):
                    rv[f'detfunnel_{sname}_{outcome}_{target}_{count}'] = 0.0
                for outcome in ('det0', 'det1', 'det2', 'det3', 'det4', 'detV'):
                    rv[f'detfunnel_{sname}_{outcome}_{target}_{count}'] = 0.0
        # record all the above keys for later use -- needs to be an explicit list
        rv['_detfunnel_keys'] = list(rv.keys())
        # Do not error-out if sim load did not work, resulting in sim_info being not present or None
        if not getattr(self, 'sim_info', False):
            return rv
        binner = RpLBins() # for is_earthlike()
        # earthlike = binner.is_earthlike(self.spc, np.arange(self.spc['nPlans']), self.spc['plan2star'])
        # hzone = binner.is_hab_zone(self.spc, np.arange(self.spc['nPlans']), self.spc['plan2star'])

        # Up here so it's available throughout
        all_stars = np.arange(self.Nstar)
        # 'promoted_stars' may not always be in the SPC dict
        promoted_stars = self.spc.get('promoted_stars', [])
        # deep-dive stars
        top_HIPs = self.sim_info['top_HIPs'] # list of hipparcos names of deep-dive stars
        top_sInds = np.where(np.isin(self.spc['Name'], top_HIPs))[0]

        ## 1: Find det status indications at each star, for each planet,
        # summing up over the DRM
        # - Status is stored as vectors-of-vectors (VoV), and represent
        #   *per-planet* counts for that star.  The VoV is a numpy (n_star,)
        #   vector of numpy "Objects", each of which is a per-planet status count.
        # - When a new det (1xNplan) is made at a star, we:
        #     ctr[sind] += this_det_status_indicator
        #   and the per-planet indicator counts will be updated.  Starting 
        #   at (scalar) zero means we don't special-case the first visit.
        # - We also need some per-star counters to keep track of visits -
        #   because planets can be [], this implies a separate counter
        n_star = self.Nstar
        # per-star counters (ordinary integers)
        tries_ctr   = np.zeros(n_star, 'int')
        fails_ctr   = np.zeros(n_star, 'int')
        success_ctr = np.zeros(n_star, 'int')
        # per-star-per-planet counters
        tries_det_ctr   = np.zeros(n_star, 'O') # VoV
        success_det_ctr = np.zeros(n_star, 'O') # VoV
        fails_det_ctr   = np.zeros(n_star, 'O') # VoV
        f_snr_det_ctr   = np.zeros(n_star, 'O') # VoV
        f_iwa_det_ctr   = np.zeros(n_star, 'O') # VoV
        f_owa_det_ctr   = np.zeros(n_star, 'O') # VoV
        for obs in self.drm:
            # skip non-detection observations
            if 'det_status' not in obs:
                continue
            sind = obs['star_ind']
            # length-nplan vector
            det_status = np.array(obs['det_status'])
            # update star counters (all are int)
            # for star:
            #     any planet succeeds => success=True
            #     no planets => success=False
            #     fail = not(success)
            success = np.any(det_status == 1)
            tries_ctr  [sind] += 1
            fails_ctr  [sind] += int(not success)
            success_ctr[sind] += int(success)
            # update planet counters (all are VoV)
            # every planet in the sind (if any) gets +1 for tries
            tries_det_ctr  [sind] += np.full(len(det_status), 1)
            success_det_ctr[sind] += (det_status ==  1)
            fails_det_ctr  [sind] += (det_status !=  1)
            f_snr_det_ctr  [sind] += (det_status ==  0)
            f_iwa_det_ctr  [sind] += (det_status == -1)
            f_owa_det_ctr  [sind] += (det_status == -2)

        # helper function, uses I[] which is modified in the loop
        I = dict()
        def select_ctr(target, star_ctr, det_ctr):
            r'''Select the counter to use depending on target.'''
            if target == 'star':
                # a scalar
                return star_ctr
            else:
                # a vector; [] if nplan = 0, no earths, etc.
                return det_ctr[I[target]]

        ## 2: Use above yields to compile metrics for promoted stars
        # can add to the "for": (top_sInds, 'deep')
        for star_list, sname in ((all_stars, 'allstar'), (promoted_stars, 'promo'), ):
            for sind in star_list:
                # quick out if not visited (all_stars case)
                if tries_ctr[sind] == 0:
                    continue
                # planets for star #sind
                plan_inds = np.where(self.spc['plan2star'] == sind)[0]
                # indicators for these planets
                I['star']    = np.full(1, True) # little used, see select_ctr
                I['allplan'] = np.full(len(plan_inds), True)
                I['hzone']   = binner.is_hab_zone (self.spc, plan_inds, sind)
                I['earth']   = binner.is_earthlike(self.spc, plan_inds, sind)

                ## 2a: candidates
                # count the candidate pool size
                #  (NB: the counter bumped only if there was a visit)
                for target in Targets:
                    rv[f'detfunnel_{sname}_cand_{target}'] += np.sum(I[target])

                ## 2b: det attempts
                # attempts are fundamentally per-star
                for target in Targets:
                    # tries is:
                    #   for star -> scalar #visits to star
                    #   for planet -> vector (1xNtargets) of visits to that planet type
                    tries = select_ctr(target, tries_ctr[sind], tries_det_ctr[sind])
                    rv[f'detfunnel_{sname}_tries_{target}_uniq'] += np.sum(tries > 0)
                    rv[f'detfunnel_{sname}_tries_{target}_cume'] += np.sum(tries)

                ## 2c: det fails
                # basically per-planet
                #    for star: all-planet-fail <=> star-fail, and no planets => fail=True
                for target in Targets:
                    # fails is just like tries above
                    fails = select_ctr(target, fails_ctr[sind], fails_det_ctr[sind])
                    # but, subcategories of "det fail" don't make sense for star
                    # (so, f_snr, etc., for star will be NaN)
                    f_snr = select_ctr(target, np.nan, f_snr_det_ctr[sind])
                    f_iwa = select_ctr(target, np.nan, f_iwa_det_ctr[sind])
                    f_owa = select_ctr(target, np.nan, f_owa_det_ctr[sind])
                    rv[f'detfunnel_{sname}_fails_{target}_uniq'] += np.sum(fails > 0)
                    rv[f'detfunnel_{sname}_fails_{target}_cume'] += np.sum(fails)
                    rv[f'detfunnel_{sname}_f_snr_{target}_cume'] += np.sum(f_snr)
                    rv[f'detfunnel_{sname}_f_iwa_{target}_cume'] += np.sum(f_iwa)
                    rv[f'detfunnel_{sname}_f_owa_{target}_cume'] += np.sum(f_owa)

                ## 2d: det success
                # basically per-planet
                #   for star: success = not(fail).
                #     any planet succeeds => success=True
                #     no planets => success=False
                for target in Targets:
                    # success is like tries above
                    # det0, det1, etc., for star still make sense
                    success = select_ctr(target, success_ctr[sind], success_det_ctr[sind])
                    rv[f'detfunnel_{sname}_success_{target}_uniq'] += np.sum(success > 0)
                    rv[f'detfunnel_{sname}_success_{target}_cume'] += np.sum(success)
                    rv[f'detfunnel_{sname}_det0_{target}_cume'] += np.sum(success == 0)
                    rv[f'detfunnel_{sname}_det1_{target}_cume'] += np.sum(success == 1)
                    rv[f'detfunnel_{sname}_det2_{target}_cume'] += np.sum(success == 2)
                    rv[f'detfunnel_{sname}_det3_{target}_cume'] += np.sum(success == 3)
                    rv[f'detfunnel_{sname}_det4_{target}_cume'] += np.sum(success == 4)
                    rv[f'detfunnel_{sname}_detV_{target}_cume'] += np.sum(success >= 5)

        ## 3: return the dict of results for this DRM
        # Note: valid keys are in '_detfunnel_keys'
        return rv


    def promotion_analysis(self):
        r'''Extract time-binned target promotion candidate numbers from DRM.

        Planet-counts are binned by instrument observing time used.
        At selected times, densities of planet-counts are also found.'''
        # Do not error-out if:
        #  (a) given an old SPC without an MsTrue field;
        #  (b) sim load did not work, resulting in sim_info being None
        # This will result in the promotion analysis not being done.
        if 'MsTrue' not in self.spc or not self.sim_info:
            return {}

        # Find "earthlike" and "hzone" for each planet in the simulation, for general use later
        binner = RpLBins()
        earthlike = binner.is_earthlike(self.spc, np.arange(self.spc['nPlans']), self.spc['plan2star'])
        hzone     = binner.is_hab_zone( self.spc, np.arange(self.spc['nPlans']), self.spc['plan2star'])
        # Determine planet period
        # M_star = stellar mass corresp. to each planet
        #   Old:
        #     M_star = 1.0*u.M_sun
        # expand MsTrue via indexing to correspond 1:1 with planets
        M_star = self.spc['MsTrue'][self.spc['plan2star']]
        # T = planet period in days
        T = np.sqrt((self.spc['a']**3 * 4 * np.pi)/(const.G * (M_star + self.spc['Mp']))).to('d')
        # T_hz{0,1} = period of {closest,farthest} HZ planet [days] (0.95,1.67 AU, rescaled by luminosity)
        # T_hz{0,1} indexed by star, because it refers to a fictional planet
        a_hz0 = 0.95*u.au * np.sqrt(self.spc['L'])
        a_hz1 = 1.67*u.au * np.sqrt(self.spc['L'])
        T_hz0 = np.sqrt((a_hz0**3 * 4 * np.pi)/(const.G * (self.spc['MsTrue'] + 0.0))).to('d')
        T_hz1 = np.sqrt((a_hz1**3 * 4 * np.pi)/(const.G * (self.spc['MsTrue'] + 0.0))).to('d')
        if 'promo-1' in DEBUG:
            print('T_hz(inner) is', T_hz0)

        # Algorithm:
        # We track 4 metrics, per-planet, vs. time:
        #   M1, count: #detections >= 3?
        #   M2, span:  detection time-span > 0.5 * planet_period
        #   M3, hzone: unique detection in the habitable zone
        #   M4, earth: unique detection of earthlike planet,
        # plus these "intersection" metrics:
        #   M3u, promo_hzone: unique planets satisfying M1 and M2 and M3.
        #   M4u, promo_earth: unique planets satisfying M1 and M2 and M4.
        # We also track a time measure, inst_time = detection_time + overhead_time,
        # which is different from wall-clock time.
        # Our output is a temporal graph, binned in time, of the number of (unique) planets
        # satisfying the above metrics at that time.

        # maintain a running tally of instrument time - detection mode only
        # these side variables are from the script, not the DRM/SPC
        #ohTime = 0.2 # [days] -- overhead time
        ohTime = self.sim_info['ohTime'] # [days] -- overhead time
        settlingTime = self.sim_info['settlingTime'] # [days] -- settling time
        my_inst_time = 0.0 # [days]
        
        # observation information: star, planets, arrival time, inst_time
        ObsInfo = namedtuple('ObsInfo', ['sind', 'pinds', 'arrival_time', 'inst_time'])
        ### [1] Accumulate per-planet and per-star detection sequences
        # planet (or star) observation times, indexed by planet number, containing list of ObsInfo tuples
        #   p_times[p] = [ObsInfo1, ..., ObsInfoN] or []
        p_times = defaultdict(list)
        s_times = defaultdict(list)
        # accumulate planet observations into p_times, star obs. into s_times
        for obs in self.drm:
            # DRM-FMT
            arrival_time = strip_units(obs['arrival_time'])
            sind = obs['star_ind']
            if 'det_time' not in obs and 'det_info' not in obs:
                continue
            # for coronagraph-only/Luvoir, detection info is kept in the 'det_info' list,
            # for starshade, detection info is in the drm entry itself (obs).
            # this abstracts the two cases by setting up a "pointer", obs_det.
            # but note, plan_inds and star_ind are always kept in obs itself.
            if 'det_info' in obs:
                obs_det = obs['det_info'][0]
            else:
                obs_det = obs
            # track the accumulated observational time
            ## NB: exoplanetObsTime includes detection + char time, we want det only, so don't use
            ## elapsed_obs_time = strip_units(obs['exoplanetObsTime'])
            my_inst_time += strip_units(obs['det_time']) + ohTime + settlingTime
            elapsed_obs_time = my_inst_time
            # detected planets
            p_detected = np.array(obs['plan_inds'])[np.where(obs_det['det_status']==1)[0]]
            # track the times we saw a planet at the star, and it could have been earthlike
            # NB: do not enforce that a *detected* planet was indeed earthlike
            if np.any(p_detected) and np.any(earthlike[obs['plan_inds']]):
                n_earth = np.sum(earthlike[obs['plan_inds']])
                s_times[sind].extend([ObsInfo(sind, obs['plan_inds'], arrival_time, elapsed_obs_time)]*n_earth)
            # track the times we saw each planet - repeats included
            # NB: a tuple is pushed: arrival_time, and cumulative-observation-time
            for p in p_detected:
                p_times[p].append(ObsInfo(sind, [p], arrival_time, elapsed_obs_time))

        if 'promo-2' in DEBUG:
            DB_p_quant = np.zeros((self.spc['nPlans'],), dtype=int)
            for p in range(self.spc['nPlans']):
                DB_p_quant[p] = binner.quantize(self.spc, p, self.spc['plan2star'][p])
            # bin boundaries -0.5, 0.5, ..., 14.5, 15.5
            DB_h_p_quant = np.histogram(DB_p_quant, bins=(np.arange(17)-0.5))[0]
            DB_earth = sum(earthlike[p] for p in p_times.keys())
            print('Total planets = %d, earths = %d, earths/stars = %.3f' % (
                earthlike.shape[0], np.sum(earthlike), np.sum(earthlike)/(1.0*self.Nstar)))
            print('Number of distinct planets detected = %d, distinct earths = %d' % (len(p_times), DB_earth))
            print('Hist:', DB_h_p_quant/np.sum(DB_h_p_quant))

        ### [2A] Compute metrics for planet promotion using p_times
        obs_thresh = 3 # this many observations is enough for the "count" metric
        # lists of times that any planet satisfied the above metrics
        # (we store the ObsInfo to allow keeping track of instrument and mission time)
        p_count_allplan, p_count_hzone, p_count_earth = [], [], []
        p_span_allplan,  p_span_hzone,  p_span_earth  = [], [], []
        p_promo_allplan, p_promo_hzone, p_promo_earth = [], [], []
        for p, p_time in six.iteritems(p_times):
            # accumulate both HZ vs. earthlike
            EARTH = earthlike[p]
            IN_HZONE = hzone[p]
            if 'promo-2' in DEBUG:
                if EARTH:
                    print('#obs of earth %d = %d' % (p, len(p_time)))
            # set up p_count --
            #   if >= obs_thresh visits: plunk on the time of the obs_thresh'th visit
            if len(p_time) >= obs_thresh:
                p_count_allplan.append(p_time[obs_thresh-1])
                if IN_HZONE:
                    p_count_hzone.append(p_time[obs_thresh-1])
                if EARTH:
                    p_count_earth.append(p_time[obs_thresh-1])
            # set up p_span, p_promo --
            #   if span > T/2: plunk on the first time it happened
            if len(p_time) > 1:
                t0 = p_time[0].arrival_time # initial arrival
                found_span_already = False
                for obs_inx, obs_info in enumerate(p_time):
                    # is the (wall clock) time-span between the first and present obs > T/2?
                    spanned = (obs_info.arrival_time - t0) > (T[p].value*0.5)
                    if spanned and not found_span_already:
                        p_span_allplan.append(obs_info)
                        if IN_HZONE:
                            p_span_hzone.append(obs_info)
                        if EARTH:
                            p_span_earth.append(obs_info)
                        found_span_already = True
                    # is the time-span sufficient, *and* >= obs_thresh visits?
                    if spanned and obs_inx >= (obs_thresh-1):
                        p_promo_allplan.append(obs_info)
                        if IN_HZONE:
                            p_promo_hzone.append(obs_info)
                        if EARTH:
                            p_promo_earth.append(obs_info)
                        break # stop now to ensure no duplicates
                        
        ### [2B] Compute per-star metrics for promotion using s_times
        # lists of times that any star satisfied the above metrics
        p_count_star, p_span_star, p_promo_star = [], [], []
        # Note: this case finds several related span metrics instead of the default, namely:
        # list of times a star satisfied several "span" metrics (depending on the T used):
        #   spanPlan:  period = min(period of planets around star) [current default]
        #   spanHZ0:   period = minimal period of a HZ planet (0.95AU)
        #   spanHZ1:   period = maximal period of a HZ planet (1.67AU)
        #   spanEarth: period = minimal period of an as-realized Earthlike planet
        p_spanPlan_star, p_spanHZ0_star, p_spanHZ1_star, p_spanEarth_star = [], [], [], []
        for s, s_time in six.iteritems(s_times):
            # set up p_count_star --
            #   if >= obs_thresh visits: plunk on the time of the obs_thresh'th visit
            if len(s_time) >= obs_thresh:
                p_count_star.append(s_time[obs_thresh-1])
            # set up p_span, p_promo --
            #   if span > T/2: plunk on the time it first happened
            if len(s_time) > 1:
                t0 = s_time[0].arrival_time # initial arrival
                # did we find a span of the above types,
                # or the default span-type ('default'), or the full promotion ('promo')
                found_span_already = {}
                for obs_inx, obs_info in enumerate(s_time):
                    # elapsed time from t0 to current obs (scalar)
                    delta_t = obs_info.arrival_time - t0
                    # 7/2019: TESTING THESE CRITERIA
                    # is the (wall clock) time-span between the first and present obs > T/2?
                    # spanPlan: any period in the star system (in effect, the smallest)
                    T_spanPlan = T[obs_info.pinds]
                    # spanHZ{0,1}: shortest/longest HZ period around star indexed s
                    T_spanHZ0 = T_hz0[s]
                    T_spanHZ1 = T_hz1[s]
                    # spanEarth: earthlike-planet periods
                    T_spanEarth = T[np.where(earthlike[obs_info.pinds])[0]]
                    # multi-pronged span check - somewhat duplicative
                    if np.any(delta_t > (0.5 * T_spanPlan.value)) and 'Plan' not in found_span_already:
                        p_spanPlan_star.append(obs_info)
                        found_span_already['Plan'] = True
                    if np.any(delta_t > (0.5 * T_spanHZ0.value)) and 'HZ0' not in found_span_already:
                        p_spanHZ0_star.append(obs_info)
                        found_span_already['HZ0'] = True
                    if np.any(delta_t > (0.5 * T_spanHZ1.value)) and 'HZ1' not in found_span_already:
                        p_spanHZ1_star.append(obs_info)
                        found_span_already['HZ1'] = True
                    if np.any(delta_t > (0.5 * T_spanEarth.value)) and 'Earth' not in found_span_already:
                        p_spanEarth_star.append(obs_info)
                        found_span_already['Earth'] = True
                    # by default: use spanPlan, the most lenient
                    spanned = np.any(delta_t > (0.5 * T_spanPlan.value))
                    if spanned and 'default' not in found_span_already:
                        p_span_star.append(obs_info)
                        found_span_already['default'] = True
                    # is the time-span sufficient, *and* >= obs_thresh visits?
                    if spanned and obs_inx >= (obs_thresh-1) and 'promo' not in found_span_already:
                        p_promo_star.append(obs_info)
                        found_span_already['promo'] = True # ensures no duplicates

        # Want to find the probability mass function of planet-count metrics at certain times
        # Compute by storing (for one sim) an indicator (0/1) vector on 0...max#planets.
        # Below: Count the first 0..cutoff-1 histogram bins, and turn it into an indicator
        # vector (with exactly 1 nonzero entry) on 0...PROMOTION_PHIST_NBINS-1
        def count_to_indicator(count, cutoff):
            return np.histogram(np.sum(count[:cutoff]), bins=PROMOTION_PHIST_BINS)[0]

        # [3] Find temporal histograms from these lists-of-times
        # Recall: p_{count,span,promo}_{allplan,hzone,earth,star} are lists of times.
        # Using these, we now, in bulk fashion:
        #   [a] bin these lists-of-times into histograms of counts-per-month;
        #   [b] make probability mass functions ("promo_phist") of these planet-counts at specific times
        # This is the list of all such histograms
        all_keys = ['%s_%s' % (k1, k2)
                        for k1 in ('count', 'span', 'promo') for k2 in ('allplan', 'hzone', 'earth', 'star')]
        # add in our family of span metrics, for the per-star case only
        all_keys.extend(['span%s_star' % metric for metric in ['Plan', 'HZ0', 'HZ1', 'Earth']])
        # compose two dictionaries holding result types [a] and [b] (maintain key order)
        promo_counts = OrderedDict()
        promo_phists = OrderedDict()
        sources = locals()
        for key in all_keys:
            # find histogram of counts, binned by time -- both instrument time, and mission clock time
            key_source = 'p_%s' % key
            hist_vs_itime = np.histogram([i.inst_time    for i in sources[key_source]], bins=PROMOTION_TIME_BINS)[0]
            hist_vs_mtime = np.histogram([i.arrival_time for i in sources[key_source]], bins=PROMOTION_TIME_BINS)[0]
            # save this "incremental" histogram and its cumulative version - instrument time
            promo_counts['h_promo_%s_incr' % key] = hist_vs_itime
            promo_counts['h_promo_%s_cume' % key] = np.cumsum(hist_vs_itime)
            if key != 'allplan':
                # converts the number of detections in "hist_vs_time" to an indicator vector over planet-counts
                # do not make this for allplanets (irrelevant, and will fall outside indicator max range)
                # one source, two destinations for two different temporal ranges (0..T_MAX)
                promo_phists['h_phist_t1_%s' % key] = count_to_indicator(hist_vs_mtime, PROMOTION_PHIST_T1_INX)
                promo_phists['h_phist_t2_%s' % key] = count_to_indicator(hist_vs_mtime, PROMOTION_PHIST_T2_INX)
                
        if 'promo-2' in DEBUG:
            print('Earth promo:')
            print(p_promo_earth)

        # [4] Compose return value
        #  -- both of the just-constructed dictionaries
        #  -- "magic" key names for each flavor of result, for use in making reductions + output
        rv = dict(_promo_count_keys = list(promo_counts.keys()),
                  _promo_phist_keys = list(promo_phists.keys()))
        rv.update(promo_counts)
        rv.update(promo_phists)
        return rv

    def per_star_promotion(self):
        r'''Compute per-star summaries relating to promotions.

        Note: All of these summaries are itemized by star.
        '''

        # for is_earthlike()
        binner = RpLBins()
        # Accumulate counts into these variables - all indexed by star number
        n_star = self.Nstar
        # count promoted stars
        promo_allplan = np.zeros(n_star)
        promo_hzone   = np.zeros(n_star)
        promo_earth   = np.zeros(n_star)

        # may not always be in the SPC dict
        promoted_stars = self.spc.get('promoted_stars', [])
        for sind in promoted_stars:
            # planets for star #sind
            plan_inds = np.where(self.spc['plan2star'] == sind)[0]
            hzones = binner.is_hab_zone (self.spc, plan_inds, sind)
            earths = binner.is_earthlike(self.spc, plan_inds, sind)
            promo_allplan[sind] = 1.0
            promo_hzone[sind]   = 1.0 * np.any(hzones)
            promo_earth[sind]   = 1.0 * np.any(earths)

        # return a dict of results for this DRM
        # Note: enforce consistency with naming conventions from here forward
        rv = {
            'h_star_promo_allplan':  promo_allplan,
            'h_star_promo_hzone':    promo_hzone,
            'h_star_promo_earth':    promo_earth,
            }
        rv['_per_star_promotion_keys'] = list(rv.keys())
        return rv


    def yield_analysis(self):
        r'''Extracts yield information from a DRM structure, for a SINGLE Exosims run.
        Uses the associated star-planet config ("spc") to classify into radius/luminosity bins
        Although the intent is to extract only yield information, some time info is found too.'''
        global VERBOSITY
        binner = RpLBins()

        # yield accumulator - a dictionary of (mostly) lists, and counters and sets, one for each band
        yac = YieldAccumulator()
        yac_bands_seen = set()

        # these are defined as sets, but they could have been Nplanet-length vectors
        set_dets_uniq    = set()
        set_chars_uniq   = set()
        set_chars_strict = set()
        # primary and alternate bin list for unique detections
        RpL_det_main = []
        RpL_det_alt  = []
        # primary and alternate bin list for all detections (mnemonic: extra detections)
        RpL_xdet_main = []
        RpL_xdet_alt  = []
        # exo-Earth counters, same scheme as RpL histograms
        exoE_det_main = 0
        exoE_det_alt = 0
        exoE_xdet_main = 0
        exoE_xdet_alt = 0
        exoE_char_strict = 0 # no SNR for strict
        exoE_char_full = 0
        exoE_char_part = 0
        exoE_xchar_full = 0 # xchar -> counts repeat chars (like xdet)
        exoE_xchar_part = 0
        SNR_exoE_char_full = []
        SNR_exoE_char_part = []
        SNR_exoE_xchar_full = []
        SNR_exoE_xchar_part = []
        # radius/luminosity for characterized targets, binned
        RpL_char_strict = []
        RpL_char_full   = []
        RpL_char_part   = []
        RpL_xchar_full  = []
        RpL_xchar_part  = []
        SNR_char_full   = []
        SNR_char_part   = []
        SNR_xchar_full   = []
        SNR_xchar_part   = []
        # list of times where detections were made
        # TODO: move to event_analysis(), so that this function only does yield
        det_time_all = []
        det_time_unq = []
        det_time_rev = []

        for obs_num, obs in enumerate(self.drm):
            plan_inds = np.array(obs['plan_inds'])
            # Process a detection:
            #   condition: 'det_time' for starshade DRMs, 'det_info' for coronagraph-only
            # DRM-FMT
            if 'det_time' in obs or 'det_info' in obs:
                # for coronagraph-only/Luvoir, detection info is kept in the 'det_info' list,
                # for starshade, detection info is in the drm entry itself (obs).
                # this abstracts the two cases by setting up a "pointer", obs_det.
                # but note, plan_inds and star_ind are always kept in obs itself.
                if 'det_info' in obs:
                    obs_det = obs['det_info'][0]
                else:
                    obs_det = obs
                det_status = obs_det['det_status']
                detections = np.where(np.array(det_status) == 1)[0]
                detected = plan_inds[detections]
                # new planet IDs at this obs.
                dets_new = set(detected).difference(set_dets_uniq)
                set_dets_uniq.update(set(detected)) # accumulate all known planet IDs
                # record detection times [day]
                arrival_time = strip_units(obs['arrival_time'])
                # all detections made
                det_time_all.extend([arrival_time] * len(detected))
                # new unique detections
                det_time_unq.extend([arrival_time] * len(dets_new))
                # revisits
                det_time_rev.extend([arrival_time] * (len(detected) - len(dets_new)))
                #print('Revisits:', (len(detected) - len(dets_new)), len(detected), len(dets_new))

                # record new detections
                for plan_id in dets_new:
                    # add in this unique detection, in whichever mode was used
                    if ('det_mode' not in obs_det) or ('combined' not in obs_det['det_mode']['instName']):
                        RpL_det_main.append(binner.quantize(self.spc, plan_id, obs['star_ind']))
                        if binner.is_earthlike(self.spc, plan_id, obs['star_ind']): exoE_det_main += 1
                    else:
                        RpL_det_alt.append( binner.quantize(self.spc, plan_id, obs['star_ind']))
                        if binner.is_earthlike(self.spc, plan_id, obs['star_ind']): exoE_det_alt += 1
                    ### TEMPORARY: a diagnostic
                    if VERBOSITY > 1 and binner.is_earthlike(self.spc, plan_id, obs['star_ind']):
                        if (exoE_det_main + exoE_det_alt) == 1: print('----') # new DRM
                        print('Earthlike planet %4d, star %4d, total -> (%2d, %2d) = %2d' % (
                            plan_id, obs['star_ind'], exoE_det_main, exoE_det_alt, exoE_det_main+exoE_det_alt))
                # record all detections
                for plan_id in detected:
                    # add in the detection, in whichever mode was used
                    if ('det_mode' not in obs_det) or ('combined' not in obs_det['det_mode']['instName']):
                        RpL_xdet_main.append(binner.quantize(self.spc, plan_id, obs['star_ind']))
                        if binner.is_earthlike(self.spc, plan_id, obs['star_ind']): exoE_xdet_main += 1
                    else:
                        RpL_xdet_alt.append( binner.quantize(self.spc, plan_id, obs['star_ind']))
                        if binner.is_earthlike(self.spc, plan_id, obs['star_ind']): exoE_xdet_alt += 1

            # Process a characterization
            if 'char_mode' in obs or 'char_info' in obs:
                # make "char_info" or a proxy of it ["char_info" is used in newer DRMs]
                #  char_info = [dict(char_time = X, char_status = Y) ...]
                if 'char_info' in obs:
                    char_info = obs['char_info']
                else:
                    char_info = [obs]
                # iterate across all detector bands that were used (e.g., 500nm, 750nm, ...)
                charizations_strict = np.array(0)
                for char in char_info:
                    char_status = char['char_status']
                    char_SNR = char['char_SNR']
                    # map from planet-id -> characterization SNR
                    plan_SNR = {p:char_SNR[i] for (i,p) in enumerate(plan_inds)}
                    # 1: The reported variable:
                    #     chars_earth_unique = exoE_char_full + exoE_char_part
                    #    it includes both full and partial chars
                    # 2: exoE_char_full will increment when *any* char_status in the char_info list
                    #    (i.e., any spectral band) has a "1"
                    #    remark: so, in multi-mode char, sometimes the red char will be partial
                    #    but a blue one would be full; this still counts as a "full" char
                    # 3: charizations_strict var insists that *all* bands (entries in char_info)
                    #    have char_status == "1"
                    #    exoE_char_strict counts these "strict" chars
                    #    So in the above scenario, a partial red char but successful blue char would
                    #    not +1 to exoE_char_strict
                    # 4: there is no strict_snr, b/c snr will differ across bands
                    # 5: set_chars_strict - strict is always "unique" and "full"
                    # strict takes the "full/part" distinction to the next level of strictness --
                    #   exoE_char_strict <= exoE_char_full <= chars_earth_unique
                    
                    charizations_strict = charizations_strict + (np.array(char_status) == 1)
                    charizations_full = np.where(np.array(char_status) ==  1)[0]
                    charizations_part = np.where(np.array(char_status) == -1)[0]
                    charized_full = plan_inds[charizations_full]
                    charized_part = plan_inds[charizations_part]

                    # per-band quantities:
                    # set_chars_uniq, chars_new_{full,part}, RpL_char_{full,part}, SNR_char_{full,part},
                    # exoE_char_{full,part}, SNR_exoE_char_{full,part}

                    for band in CHAR_BANDS:
                        if not char_within_band(char, band): continue
                        yac_bands_seen.add(band)
                        _band = '_' + band # for ease of naming
                        # full chars in this band
                        chars_new_full = set(charized_full).difference(yac['set_chars_uniq'+_band])
                        for plan_id in chars_new_full:
                            yac['RpL_char_full'+_band].append(binner.quantize(self.spc, plan_id, obs['star_ind']))
                            yac['SNR_char_full'+_band].append(plan_SNR[plan_id])
                            if binner.is_earthlike(self.spc, plan_id, obs['star_ind']):
                                yac['exoE_char_full'+_band] += 1
                                yac['SNR_exoE_char_full'+_band].append(plan_SNR[plan_id])
                        # partial chars in this band
                        chars_new_part = set(charized_part).difference(yac['set_chars_uniq'+_band])
                        for plan_id in chars_new_part:
                            yac['RpL_char_part'+_band].append(binner.quantize(self.spc, plan_id, obs['star_ind']))
                            yac['SNR_char_part'+_band].append(plan_SNR[plan_id])
                            if binner.is_earthlike(self.spc, plan_id, obs['star_ind']):
                                yac['exoE_char_part'+_band] += 1
                                yac['SNR_exoE_char_part'+_band].append(plan_SNR[plan_id])
                        # keep a running tabulation of all characterizations so far in this band
                        yac['set_chars_uniq'+_band].update(set(charized_full))
                        yac['set_chars_uniq'+_band].update(set(charized_part))

                    # new chars for this obs, found by removing all prior chars (full OR partial)
                    chars_new_full = set(charized_full).difference(set_chars_uniq)
                    chars_new_part = set(charized_part).difference(set_chars_uniq)
                    # binned characterizations (new at this obs), full and partial
                    # plus, the characterization SNR corresponding to each
                    for plan_id in chars_new_full:
                        RpL_char_full.append(binner.quantize(self.spc, plan_id, obs['star_ind']))
                        SNR_char_full.append(plan_SNR[plan_id])
                        if binner.is_earthlike(self.spc, plan_id, obs['star_ind']):
                            exoE_char_full += 1
                            SNR_exoE_char_full.append(plan_SNR[plan_id])
                    for plan_id in chars_new_part:
                        RpL_char_part.append(binner.quantize(self.spc, plan_id, obs['star_ind']))
                        SNR_char_part.append(plan_SNR[plan_id])
                        if binner.is_earthlike(self.spc, plan_id, obs['star_ind']):
                            exoE_char_part += 1
                            SNR_exoE_char_part.append(plan_SNR[plan_id])
                    # binned chars (new or repeat at this obs), full + partial
                    # (note, xchar is analogous to xdet)
                    for plan_id in charized_full:
                        RpL_xchar_full.append(binner.quantize(self.spc, plan_id, obs['star_ind']))
                        SNR_xchar_full.append(plan_SNR[plan_id])
                        if binner.is_earthlike(self.spc, plan_id, obs['star_ind']):
                            exoE_xchar_full += 1
                            SNR_exoE_xchar_full.append(plan_SNR[plan_id])
                    for plan_id in charized_part:
                        RpL_xchar_part.append(binner.quantize(self.spc, plan_id, obs['star_ind']))
                        SNR_xchar_part.append(plan_SNR[plan_id])
                        if binner.is_earthlike(self.spc, plan_id, obs['star_ind']):
                            exoE_xchar_part += 1
                            SNR_exoE_xchar_part.append(plan_SNR[plan_id])
                    # keep a running tabulation of all characterizations so far
                    set_chars_uniq.update(set(charized_full))
                    set_chars_uniq.update(set(charized_part))
                # (end of for-loop over detector integration bands -- handle strict chars)
                # planets that were characterized in every one of the detector bands above
                charized_strict = plan_inds[charizations_strict == len(char_info)]
                chars_new_strict = set(charized_strict).difference(set_chars_strict)
                for plan_id in chars_new_strict:
                    RpL_char_strict.append(binner.quantize(self.spc, plan_id, obs['star_ind']))
                    if binner.is_earthlike(self.spc, plan_id, obs['star_ind']):
                        exoE_char_strict += 1 # (no SNR in this case)
                set_chars_strict.update(set(charized_strict))

        # Earth histograms -- just a single count, actually
        h_earth_char_all    = np.histogram(np.array([exoE_char_full + exoE_char_part]),
                                            EARTH_CHAR_COUNT_BINS)[0]
        h_earth_xchar_all   = np.histogram(np.array([exoE_xchar_full + exoE_xchar_part]),
                                            EARTH_CHAR_COUNT_BINS)[0]
        h_earth_char_strict = np.histogram(np.array([exoE_char_strict]),
                                            EARTH_CHAR_COUNT_BINS)[0]

        # Find radius/luminosity histograms
        RpL_bin_edges = binner.RpL_bin_edge_list
        # [1] unique detections
        h_RpL_det_main = np.histogram(RpL_det_main, RpL_bin_edges)[0]
        h_RpL_det_alt  = np.histogram(RpL_det_alt,  RpL_bin_edges)[0]
        # [2] all detections
        h_RpL_xdet_main = np.histogram(RpL_xdet_main, RpL_bin_edges)[0]
        h_RpL_xdet_alt  = np.histogram(RpL_xdet_alt,  RpL_bin_edges)[0]
        # [3] full-population histograms -- as a check on the simulated universe/target list
        # [3a] for the RpL bins
        nPlans = self.spc['nPlans']
        # quantizer is not vectorized, so must loop
        RpL_population = np.zeros((nPlans,), dtype=int)
        for p in range(nPlans):
            RpL_population[p] = binner.quantize(self.spc, p, self.spc['plan2star'][p])
        # planet counts within each bin / #stars
        # NB: there can be out-of-range planets in the Rp,L space
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning) # ignore 0/0 -> nan
            h_RpL_population = np.histogram(RpL_population, RpL_bin_edges)[0] / (1.0 * self.Nstar)
        # [3b] for Earthlike; this the empirical eta-Earth
        earthlike = binner.is_earthlike(self.spc, np.arange(nPlans), self.spc['plan2star'])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning) # ignore 0/0 -> nan
            exoE_population = np.sum(earthlike) / (1.0 * self.Nstar)

        # returned value
        rv = {}
        # per-band returned values
        for band in yac_bands_seen:
            _band = '_' + band
            # radius/luminosity histograms: characterization
            rv['h_RpL_char_full'+_band] = np.histogram(yac['RpL_char_full'+_band], RpL_bin_edges)[0]
            rv['h_RpL_char_part'+_band] = np.histogram(yac['RpL_char_part'+_band], RpL_bin_edges)[0]
            # characterization SNR
            rv['h_RpL_char_snr' +_band] = stats.binned_statistic(x=yac['RpL_char_full'+_band],
                                                    values=yac['SNR_char_full'+_band],
                                                    statistic='mean',
                                                    bins=RpL_bin_edges,
                                                    range=(RpL_bin_edges[0],RpL_bin_edges[-1]))[0]
            # exo-Earth counts
            rv['exoE_char_full'+_band] = yac['exoE_char_full'+_band]
            rv['exoE_char_part'+_band] = yac['exoE_char_part'+_band]
            # exo-Earth SNR
            snr_list = yac['SNR_exoE_char_full'+_band]
            rv['exoE_char_snr' +_band] = np.mean(snr_list) if snr_list else np.nan

        # radius/luminosity histograms: characterization
        #  this histogram drops the "0", or out-of-range, RpL region
        h_RpL_char_strict = np.histogram(RpL_char_strict, RpL_bin_edges)[0]
        h_RpL_char_full   = np.histogram(RpL_char_full,   RpL_bin_edges)[0]
        h_RpL_char_part   = np.histogram(RpL_char_part,   RpL_bin_edges)[0]
        h_RpL_xchar_full  = np.histogram(RpL_xchar_full,  RpL_bin_edges)[0]
        h_RpL_xchar_part  = np.histogram(RpL_xchar_part,  RpL_bin_edges)[0]
        # throughput: #(char) / #(there)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning) # ignore 0/0 -> nan
            h_RpL_char_tput_strict = h_RpL_char_strict / (1.0*self.Nstar*h_RpL_population)
            h_RpL_char_tput_full   = h_RpL_char_full   / (1.0*self.Nstar*h_RpL_population)
            h_RpL_xchar_tput_full  = h_RpL_xchar_full  / (1.0*self.Nstar*h_RpL_population)
            # 10/2020: 0/0 has been promoted from RuntimeWarning to ZeroDivisionError, must special-case
            if np.sum(earthlike) > 0:
                exoE_char_tput_strict  = exoE_char_strict  / (1.0 * np.sum(earthlike)) # scalar, nan OK
                exoE_char_tput_full    = exoE_char_full    / (1.0 * np.sum(earthlike)) # scalar, nan OK
                exoE_xchar_tput_full   = exoE_xchar_full   / (1.0 * np.sum(earthlike)) # scalar, nan OK
            else:
                exoE_char_tput_strict  = np.nan * exoE_char_strict
                exoE_char_tput_full    = np.nan * exoE_char_full
                exoE_xchar_tput_full   = np.nan * exoE_xchar_full
        # compute the mean of SNRs within each of the RpL characterization bins
        #   Note: empty bins are filled with NaN
        #   Note: this is using only full characterizations
        #   Note: the explicit range cutoff seems unnecessary, but the python docs are ambiguous
        #   Note: Later on, we find the mean (across DRMs) of *these* means.
        h_RpL_char_snr = stats.binned_statistic(x=RpL_char_full,
                                                    values=SNR_char_full,
                                                    statistic='mean',
                                                    bins=RpL_bin_edges,
                                                    range=(RpL_bin_edges[0],RpL_bin_edges[-1]))[0]
        h_RpL_xchar_snr = stats.binned_statistic(x=RpL_xchar_full,
                                                    values=SNR_xchar_full,
                                                    statistic='mean',
                                                    bins=RpL_bin_edges,
                                                    range=(RpL_bin_edges[0],RpL_bin_edges[-1]))[0]
        # scalar for exo-earths
        exoE_char_snr  = np.mean(SNR_exoE_char_full)  if SNR_exoE_char_full  else np.nan
        exoE_xchar_snr = np.mean(SNR_exoE_xchar_full) if SNR_exoE_xchar_full else np.nan

        # FIXME: NO LONGER USED, REMOVE
        # bin the detection-times ("h_" is mnemonic for histogrammed)
        ## h_det_time_all = np.histogram(det_time_all, DETECTION_TIME_BINS)[0]
        ## h_det_time_unq = np.histogram(det_time_unq, DETECTION_TIME_BINS)[0]
        ## h_det_time_rev = np.histogram(det_time_rev, DETECTION_TIME_BINS)[0]

        # some portions of the return value are automated
        namespace = locals()
        # Rp/L histograms -- "radlum" family
        qoi_radlum = [
            'h_RpL_det_main',
            'h_RpL_det_alt',
            'h_RpL_xdet_main',
            'h_RpL_xdet_alt',
            'h_RpL_char_strict',
            'h_RpL_char_full',
            'h_RpL_char_part',
            'h_RpL_char_snr',
            'h_RpL_xchar_full',
            'h_RpL_xchar_part',
            'h_RpL_xchar_snr',
            'h_RpL_population',
            'h_RpL_char_tput_strict',
            'h_RpL_char_tput_full',
            'h_RpL_xchar_tput_full',
            ]
        rv_radlum = {qoi: namespace[qoi] for qoi in qoi_radlum}
        rv.update(rv_radlum)
        rv['_radlum_keys'] = qoi_radlum
        # counts for exo-Earths -- "earth" family
        qoi_earth = [
            'exoE_det_main',
            'exoE_det_alt',
            'exoE_xdet_main',
            'exoE_xdet_alt',
            'exoE_char_strict',
            'exoE_char_full',
            'exoE_char_part',
            'exoE_xchar_full',
            'exoE_xchar_part',
            'exoE_char_snr',
            'exoE_xchar_snr',
            'exoE_population',
            'exoE_char_tput_strict',
            'exoE_char_tput_full',
            'exoE_xchar_tput_full',
            ]
        rv_earth = {qoi: namespace[qoi] for qoi in qoi_earth}
        rv.update(rv_earth)
        rv['_earth_keys'] = qoi_earth
        # a bag of yield figures - reduction is NOT handled the way as other QOIs
        rv0 = {
            'dets_unique':     np.array(list(set_dets_uniq)),
            'chars_unique':    np.array(list(set_chars_uniq)),
            'chars_strict':    np.array(list(set_chars_strict)),
            # which bands (red/blue) were seen in this DRM
            'char_bands_seen': yac_bands_seen, # a set
            }
        rv.update(rv0)
        # more yield figures - reduced like everything else
        rv1 = {
            # histogram for exo-Earths -- earth_char_count family
            'h_earth_char_all':    h_earth_char_all,
            'h_earth_xchar_all':   h_earth_xchar_all,
            'h_earth_char_strict': h_earth_char_strict,
            # times -- part of  "times" family, but other things there too
            ## 'h_det_time_all': h_det_time_all,
            ## 'h_det_time_unq': h_det_time_unq,
            ## 'h_det_time_rev': h_det_time_rev, 
            }
        # return the pooled result
        rv.update(rv1)
        return rv

    def yield_time_analysis(self):
        r'''Extracts yield-vs-time information from a DRM structure, for a SINGLE Exosims run.

        Uses the associated star-planet config ("spc") to identify exo-Earths.
        This method handles only mission-time-related detection/characterization information, 
        while yield_analysis() handles non-time information.'''
        global VERBOSITY
        binner = RpLBins()

        # yield accumulator - a dictionary of (mostly) lists, and counters and sets, one for each band
        yac = YieldAccumulator()
        yac_bands_seen = set()

        # these are defined as sets, but they could have been Nplanet-length vectors
        set_dets_uniq = set()
        set_chars_uniq = set()

        for obs_num, obs in enumerate(self.drm):
            plan_inds = np.array(obs['plan_inds'])
            arrival_time = strip_units(obs['arrival_time']) # [day]
            # Process a detection:
            #   condition: 'det_time' for starshade DRMs, 'det_info' for coronagraph-only
            # DRM-FMT
            if 'det_time' in obs or 'det_info' in obs:
                # for coronagraph-only/Luvoir, detection info is kept in the 'det_info' list,
                # for starshade, detection info is in the drm entry itself (obs).
                # this abstracts the two cases by setting up a "pointer", obs_det.
                # but note, plan_inds and star_ind are always kept in obs itself.
                if 'det_info' in obs:
                    obs_det = obs['det_info'][0]
                else:
                    obs_det = obs
                det_status = obs_det['det_status']
                detections = np.where(np.array(det_status) == 1)[0]
                detected = plan_inds[detections]
                # new planet IDs at this obs.
                dets_new = set(detected).difference(set_dets_uniq)
                set_dets_uniq.update(set(detected)) # accumulate all known planet IDs
                # record cumulative, unique, and revisit detections, for all-planets and Earths
                for plan_id in detected:
                    is_earth = binner.is_earthlike(self.spc, plan_id, obs['star_ind'])
                    yac['time_det_allplan_cume'].append(arrival_time)
                    if is_earth:
                        yac['time_det_earth_cume'].append(arrival_time)
                    if plan_id in dets_new:
                        # add in this unique detection
                        yac['time_det_allplan_uniq'].append(arrival_time)
                        if is_earth:
                            yac['time_det_earth_uniq'].append(arrival_time)
                    else:
                        yac['time_det_allplan_revi'].append(arrival_time)
                        if is_earth:
                            yac['time_det_earth_revi'].append(arrival_time)

            # Process a characterization
            if 'char_mode' in obs or 'char_info' in obs:
                # make "char_info" or a proxy of it ["char_info" is used in newer DRMs]
                #  char_info = [dict(char_time = X, char_status = Y) ...]
                if 'char_info' in obs:
                    char_info = obs['char_info']
                else:
                    char_info = [obs]
                # accumulate across multiple bands (e.g., red + blue)
                for char in char_info:
                    char_status = char['char_status']
                    charizations_full = np.where(np.array(char_status) ==  1)[0]
                    charizations_part = np.where(np.array(char_status) == -1)[0]
                    charized_full = plan_inds[charizations_full]
                    charized_part = plan_inds[charizations_part]
                    # for this plot, part = strict_partial UNION full
                    charized_part = np.append(charized_part, charized_full)
                    # per-band quantities: set_chars_uniq_{band}
                    # time_char_{full,part}_{allplan,earth}_{cume,uniq,revi}_{band}
                    for band in CHAR_BANDS:
                        if not char_within_band(char, band): continue
                        yac_bands_seen.add(band)
                        _band = '_' + band # for ease of naming
                        # new full chars *in this band*
                        chars_new_full = set(charized_full).difference(yac['set_chars_full_uniq'+_band])
                        for plan_id in charized_full:
                            is_earth = binner.is_earthlike(self.spc, plan_id, obs['star_ind'])
                            yac['time_char_full_allplan_cume'+_band].append(arrival_time)
                            if is_earth:
                                yac['time_char_full_earth_cume'+_band].append(arrival_time)
                            if plan_id in chars_new_full:
                                yac['time_char_full_allplan_uniq'+_band].append(arrival_time)
                                if is_earth:
                                    yac['time_char_full_earth_uniq'+_band].append(arrival_time)
                            else:
                                yac['time_char_full_allplan_revi'+_band].append(arrival_time)
                                if is_earth:
                                    yac['time_char_full_earth_revi'+_band].append(arrival_time)
                        # partial chars *in this band*
                        # NB: if full char was done already, this is not a new partial char
                        chars_new_part = set(charized_part).difference(yac['set_chars_part_uniq'+_band])
                        for plan_id in charized_part:
                            is_earth = binner.is_earthlike(self.spc, plan_id, obs['star_ind'])
                            yac['time_char_part_allplan_cume'+_band].append(arrival_time)
                            if is_earth:
                                yac['time_char_part_earth_cume'+_band].append(arrival_time)
                            if plan_id in chars_new_part:
                                yac['time_char_part_allplan_uniq'+_band].append(arrival_time)
                                if is_earth:
                                    yac['time_char_part_earth_uniq'+_band].append(arrival_time)
                            else:
                                yac['time_char_part_allplan_revi'+_band].append(arrival_time)
                                if is_earth:
                                    yac['time_char_part_earth_revi'+_band].append(arrival_time)

                        # keep a running tabulation of all characterizations so far in this band
                        yac['set_chars_full_uniq'+_band].update(set(charized_full))
                        yac['set_chars_part_uniq'+_band].update(set(charized_part))

        # name of each time-list within yac[] to convert to a histogram
        names = []
        names.extend(['time_det_%s_%s' % (target, status)
                        for target in ('allplan', 'earth')
                        for status in ('cume', 'uniq', 'revi')])
        names.extend(['time_char_%s_%s_%s_%s' % (success, target, status, band)
                        for success in ('full', 'part')
                        for target in ('allplan', 'earth')
                        for status in ('cume', 'uniq', 'revi')
                        for band in yac_bands_seen])
                    
        # accumulate the returned value (rv) by binning the detection-times
        # ("h_" is mnemonic for histogrammed)
        rv = {}
        for n in names:
            rv['h_' + n] = np.histogram(yac[n], DETECTION_TIME_BINS)[0]
        # the list of keys that we're returning ... for use in later steps
        rv['_yield_time_keys'] = ['h_' + n for n in names]

        # make an index of target depletion using the events-vs-time accumulator
        # that we already have; they will be injected in the final result
        rv2 = self.target_depletion_subanalysis(names, yac)
        rv2['_target_depletion_keys'] = list(rv2.keys())

        # return the pooled result
        return {**rv, **rv2}


    def target_depletion_subanalysis(self, names, yac):
        r'''
        Uses the yield accumulator (list of times of yield events) to
        find metrics for target depletion
        Note: the most chacteristic field of interest is:
          'time_char_%s_%s_%s_%s' % ('full', 'allplan', 'uniq', 'union')
        '''
        rv = dict()
        for n in names:
            # 1: Filter. These keys will not be used to find any
            # target depletion value, not even a placeholder.
            # filter down to just chars
            if 'time_char' not in n:
                continue # skip detections
            # filter out some others that don't seem interesting
            # revisits, red/blue bands
            if ('_revi' in n or
                '_red' in n or
                '_blue' in n):
                continue
            # 2: Compute metrics. All keys below here will have a metric.
            # define the names: n_orig is the source name, and we strip off
            # the leading time_ to make the destination name
            n_orig = n
            n_dest = n[5:]
            # no yield at all in that category => yac[n_orig] will not have anything
            # filled in - but yac is a list-accumulator, so the slot will materialize
            # as empty. Fill in a dummy value if the materialized list is [].
            if len(yac[n_orig]) == 0:
                slope0, slope1, slope2 = np.array(0.0), np.array(0.0), np.array(0.0)
                t80 = np.array(np.nan)
            else:
                ycn = yac[n_orig]
                # TODO: key off missionLife instead
                # ratio: (last year yield) / (full yield), and friends
                c_sta_1 = len([y1 for y1 in ycn if y1 < (1*365.25)])
                c_fin_1 = len([y1 for y1 in ycn if y1 > (DETECTION_TIME_BINS[-1] - 1*365.25)])
                c_fin_2 = len([y1 for y1 in ycn if y1 > (DETECTION_TIME_BINS[-1] - 2*365.25)])
                c_total = len(ycn)
                slope0 = np.array(c_sta_1 / c_total) if c_total > 0 else np.array(0.0)
                slope1 = np.array(c_fin_1 / c_total) if c_total > 0 else np.array(0.0)
                slope2 = np.array((c_fin_2-c_fin_1) / c_total) if c_total > 0 else np.array(0.0)
                # time-to-yield: mission time where we reach 80% of full yield
                # t80: we round down, and we are not interpolating -- if yield = 7
                # then n80 = 0.80 * 7 = 5.6 -> 5 and we take the time of the
                # 5'th successful char, at index = 4.
                # If yield = 1, n80 = 0.8 -> 0 and we want index (-1). We cheat
                # this one and take the time of the first and only char.
                N_yield = len(ycn)
                n80 = int(np.floor(0.80 * N_yield))
                t80 = np.array(ycn[max(n80 - 1, 0)])
            rv[f'tdep_slope_yp1_{n_dest}'] = slope0 # year 1
            rv[f'tdep_slope_ym1_{n_dest}'] = slope1 # final year
            rv[f'tdep_slope_ym2_{n_dest}'] = slope2 # final-but-one year
            rv[f'tdep_t80_{n_dest}'] = t80 # time-to-80%
        return rv


    def visit_time_analysis(self):
        r'''Extracts visits-vs-time information from a DRM structure, for a SINGLE Exosims run.

        Method handles mission-time-related detection/characterization information,
        counting by stars (visits) not counting by planets (yield).'''
        global VERBOSITY

        # visit accumulator - a dictionary of lists of arrival-times
        names = [f'{typ}_{status}'
                     for typ in ('det', 'char')
                     for status in ('visit', 'revi', 'uniq')]
        vac = {key: [] for key in names}

        # allow us to track revisits
        stars_visited_det  = set()
        stars_visited_char = set()

        for obs_num, obs in enumerate(self.drm):
            sind = obs['star_ind']
            plan_inds = np.array(obs['plan_inds'])
            arrival_time = strip_units(obs['arrival_time']) # [day]
            # Process a detection:
            #   condition: 'det_time' for starshade DRMs, 'det_info' for coronagraph-only
            if 'det_time' in obs or 'det_info' in obs:
                # record cumulative, unique, and revisits
                vac['det_visit'].append(arrival_time)
                if sind in stars_visited_det:
                    vac['det_revi'].append(arrival_time)
                else:
                    vac['det_uniq'].append(arrival_time)
                    stars_visited_det.add(sind)
            # Process a characterization
            if 'char_mode' in obs or 'char_info' in obs:
                # record cumulative, unique, and revisits
                vac['char_visit'].append(arrival_time)
                if sind in stars_visited_char:
                    vac['char_revi'].append(arrival_time)
                else:
                    vac['char_uniq'].append(arrival_time)
                    stars_visited_char.add(sind)

        # accumulate the returned value (rv) by binning the visit-times
        # ("h_" is mnemonic for histogrammed)
        rv = {}
        for n in names:
            rv['h_visit_' + n] = np.histogram(vac[n], DETECTION_TIME_BINS)[0]
        # the list of keys that we're returning ... for use in later steps
        rv['_visit_time_keys'] = list(rv.keys())
        # return the pooled result
        return rv


    def summarize(self, econo=True):
        r'''Find the summary of the sim as a dictionary held within the object.

        The convention is that the summary is built up by calling a series of analytic
        routines, each of which returns a dictionary of summary information.  The overall
        summary dictionary is a union of each individual summary.
        If econo, delete the DRM and keep only the summary.'''
        # this dict holds reductions for the current sim
        summary = {}
        # fold in summary of yield (det/char, all/exo-earth, red/blue)
        summary.update(self.yield_analysis())
        # fold in summary of yield-vs-time (det/char, all/exo-earth, red/blue)
        summary.update(self.yield_time_analysis())
        # fold in summary of resource use (fuel, integration time)
        summary.update(self.resource_analysis())
        # fold in delta-v vs time (related to fuel)
        summary.update(self.delta_v_analysis())
        # fold in star-visit timeline
        summary.update(self.visit_time_analysis())
        # fold in summary of revisits-vs-time
        summary.update(self.summarize_revisits())
        # fold in per-star summary of yield, tInt
        summary.update(self.per_star_yield())
        # fold in event durations
        summary.update(self.event_analysis())
        # fold in event counts
        summary.update(self.count_events())
        # fold in target promotion summaries
        summary.update(self.promotion_analysis())
        # fold in per-star promotion summary
        summary.update(self.per_star_promotion())
        # fold in "funnel" promotion and deep-dive summaries
        summary.update(self.funnel_analysis())
        # fold in "detfunnel" detection summary
        summary.update(self.det_funnel_analysis())
        # delete the base data if asked
        if econo:
            self.drm = None
            self.spc = None
        # keep a reference in the object (possibly not advisable?)
        self.summary = summary
        # also return the summary-dictionary
        return summary


# NB: Present in outer scope as a helper for load_and_reduce() below.
def outer_load_and_reduce(f, args=None):
    r'''Load a sim and summarize it into a dict.

    This must be present at the outer scope of the file so it can be loaded
    by a separate process that is created by the multiprocessing module.'''
    if args is not None and args.verbose > 0:
        print('Processing <%s> in pid #%d' % (f, os.getpid()))
    sim = SimulationRun(f, sim_info=args.sim_info if args else {})
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
        self.auto_keys = dict() # automatically handled keys for each analysis phase
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

    def get_key_family_attrs(self, reductions, key):
        r'''Pool the attrs named by the given key (the key family) across all reductions.'''
        attr_set = set()
        for r in reductions:
            # r.get(): allow processing to continue if the family was not present
            # 2024: but, processing still may not continue to completion
            attr_set.update(set(r.get(key, [])))
        # ensure the order is deterministic
        return list(sorted(attr_set))

    def get_attr_null_value(self, attr, valid_value):
        r'''Utility, return the null value for the attribute named attr.

        The null value is the value which, when averaged in with valid values, will produce
        the correct ensemble mean in the case where the underlying run didn't have that
        attribute.  Not having an attribute happens currently for multi-band char, for which 
        some runs do not have any chars done in the red band (for example), but other runs
        do.  Affected attributes seen in this routine are:

           exoE_char_{full,part,snr}_B
           h_RpL_char_{full,part,snr}_B
           h_time_char_{full,part}_{allplan,earth}_{cume,revi}_B

        where, B is a band that could be 'union' or 'red'.  In this case, the
        exoE_* parameters are float's, while the h_time_... are np.ndarray's.
        '''
        if isinstance(valid_value, float):
            # float or np.float -- e.g., exoEarth char count or SNR
            if 'snr' in attr:
                return np.nan # NaN's are not counted in the mean, it's formally a 0/0
            else:
                return 0.0
        elif isinstance(valid_value, np.ndarray):
            # this will have the correct shape
            if 'snr' in attr:
                return np.nan * valid_value
            else:
                return 0.0 * valid_value # this is non-NaN
        else:
            assert False, 'Unknown type for missing value in %s' % attr

    def reduce_only(self):
        r'''Make reductions as a list-of-dicts, one per sim.  

        NOTE: This routine is obsolete.
        We now prefer to load-and-reduce each sim, rather than loading all sims 
        and then reducing them all.
        '''
        reductions = [sim.summarize() for sim in self.sims]
        # re-group, and then accumulate the above reductions across DRMs
        # result is a dict containing reduced data, stored as self.summary
        self.regroup_and_accum(reductions)
        # set up pass-thru entries so they can be dumped with the rest
        self.insert_static_entries()

    def load_and_reduce(self):
        r'''Load sim drm/spc, reduce each sim, accumulate summaries across sims.

        Each dict gives one summary statistic over that single sim.'''

        # In general, creates a pool of workers (separate unix processes)
        # (but: if jobs <= 1, it uses ordinary python map() and does no multiprocessing)
        with WorkerMap(self.args.jobs) as map_function:
            # map the load-and-reduce function over each file
            # reductions is a list of dicts containing summaries
            # Note: for python3, need to ensure we materialize the list
            reductions = list(map_function(partial(outer_load_and_reduce, args=self.args),
                                      self.sim_files))
        # hacky fix for no-drm case
        if len(self.sim_files) == 0:
            reductions = [outer_load_and_reduce(None)]
        # re-group, and then accumulate the above reductions across sims
        # result is a dict containing reduced data, stored as self.summary
        self.regroup_and_accum(reductions)
        # set up pass-thru entries so they can be dumped with the rest
        self.insert_static_entries()

    def insert_static_entries(self):
        r'''Insert extra values, that are pass-thru's from other sources, to the summary.

        The distinction from the other values is that these values are not means, 
        quantiles, et., across the DRM ensemble.'''
        # start empty
        summary = {}
        # A1: detection time bins, lo + hi
        summary['h_det_time_lo'] = DETECTION_TIME_BINS[:-1]
        summary['h_det_time_hi'] = DETECTION_TIME_BINS[1:]
        # A2: resolutions of event duration bins (same length)
        summary['h_event_b0_duration_lo'] = DURATION_TIME_B0_BINS[:-2]
        summary['h_event_b0_duration_hi'] = DURATION_TIME_B0_BINS[1:-1]
        summary['h_event_b1_duration_lo'] = DURATION_TIME_B1_BINS[:-2]
        summary['h_event_b1_duration_hi'] = DURATION_TIME_B1_BINS[1:-1]
        summary['h_event_b2_duration_lo'] = DURATION_TIME_B2_BINS[:-2]
        summary['h_event_b2_duration_hi'] = DURATION_TIME_B2_BINS[1:-1]
        # A3: event count bins
        summary['h_event_count_lo'] = np.arange(EVENT_COUNT_NBINS)
        summary['h_event_count_hi'] = np.arange(EVENT_COUNT_NBINS) + 1
        # A4: promotion histogram bins
        summary['h_promo_time_lo'] = PROMOTION_TIME_BINS[:-1]
        summary['h_promo_time_hi'] = PROMOTION_TIME_BINS[1:]
        # A5: promotion planet-histogram bins
        summary['h_phist_count_lo'] = PROMOTION_PHIST_BINS[:-1]
        summary['h_phist_count_hi'] = PROMOTION_PHIST_BINS[1:]
        # B: Rp/L histogram low/hi bin boundaries
        #     note, this parameterization is the same size as the data,
        #     so it works well with CSV
        summary['h_Rp_lo'] = RpLBins.Rp_lo.ravel()
        summary['h_Rp_hi'] = RpLBins.Rp_hi.ravel()
        summary['h_L_lo']  = RpLBins.L_lo.ravel()
        summary['h_L_hi']  = RpLBins.L_hi.ravel()
        # B2: earth char count bins
        summary['h_earth_char_count_lo'] = EARTH_CHAR_COUNT_BINS[:-1]
        summary['h_earth_char_count_hi'] = EARTH_CHAR_COUNT_BINS[1:]

        # C: visit-count bins
        #    again, keep it the same size as the data (data is found elsewhere)
        #    this is generated as a special case of the basic summarizer-function
        dummy_sim = SimulationRun(None)
        summary['h_visit_bins'] = dummy_sim.summarize_revisits()['h_visit_bins']
        # D: star characteristics
        #    one value will be supplied for each characteristic, for each star
        #    we pass each property through prefixed by "star_"
        star_attrs = [('Name', 'star_Name', np_force_string),
                      ('L',    'star_L',    strip_units),
                      ('dist', 'star_dist', strip_units),
                      ('Spec', 'star_Spec', np_force_string)]
        if self.Ndrm_actual > 0:
            # there was at least one real DRM -> at least one real SPC
            # it was saved within the ensemble
            spc1 = self.spc
            for attr_old, attr_new, converter in star_attrs:
                summary[attr_new] = converter(spc1[attr_old])
        else:
            # put something there
            for attr_old, attr_new, converter in star_attrs:
                summary[attr_new] = []
        summary['star_sInd'] = np.arange(self.Nstar)
        # place the entries in the root object
        self.summary.update(summary)

    def regroup_and_accum(self, reductions):
        r'''Accumulate various summaries across the ensemble.

        Nothing is returned: result is placed in the object state.'''
        
        # 0a: Make the list of attributes to accumulate
        # TODO: eliminate the copy-pasted key-lists such as below
        # by using the auto-keys mechanism in (0b) below
        attrs = [
            # NB: these 3 attr's are computed in yield_analysis()
            # 'h_det_time_all', 'h_det_time_unq', 'h_det_time_rev',
            # attr's computed in yield_analysis
            # earth-char histograms as a function of #chars
            'h_earth_char_all', 'h_earth_xchar_all', 'h_earth_char_strict',
            ]
        # 0b: add in auto-generated keys, using a key family name convention
        #     - this flows information from reductions into attrs and self.auto_keys
        #     - The "auto_keys" property is used later, when writing values out.
        #     - Note: not all keys may be common across all reductions, because
        #     of dynamic keys (multi-band chars).  So we pool keys across all reductions.
        #     This affects '_yield_time_keys', for instance.
        attrs_auto = []
        for key_family in ('event_count', 'event_analysis', 'promo_count', 'promo_phist',
                               'funnel', 'detfunnel', 
                               'per_star_yield', 'per_star_promotion', 'resource_analysis',
                               'summarize_revisits', 'visit_time', 
                               'earth_char', 'yield_time', 'target_depletion', 
                               'radlum', 'earth', 'delta_v'):
            # e.g., all 'promo_count' attrs are found by looking up '_promo_count_keys' in reductions
            attrs_new = self.get_key_family_attrs(reductions, '_%s_keys' % key_family)
            attrs_auto.extend(attrs_new)
            self.auto_keys[key_family] = attrs_new
        attrs.extend(attrs_auto)
            
        # FIXME: replace this with auto-generated keys, as above
        # 0c: add dynamic attrs due to per-band chars
        # determine which bands are present
        char_bands_seen = set()
        for r in reductions:
            char_bands_seen.update(r['char_bands_seen'])
        char_bands_seen = list(sorted(char_bands_seen)) # make the order deterministic
        for band in char_bands_seen:
            # per-band chars
            for attr in ('h_RpL_char_full', 'h_RpL_char_part', 'h_RpL_char_snr'):
                attrs.append(attr + '_' + band)
            # per-band chars, exo-Earth
            for attr in ('exoE_char_full', 'exoE_char_part', 'exoE_char_snr'):
                attrs.append(attr + '_' + band)
        # 1: transpose the reductions from [drm][attribute] to [attribute][drm]
        #   Then, below, we then flatten along the final index, [drm].
        #   Because of dynamic attrs above, not all attrs may be in all DRMs
        # accum: dictionary of lists, indexed by attr
        accum = defaultdict(list)
        place_holder = {'Value': 'Unknown'} # this token is inserted when the reduction is not present
        # skipped: dictionary of run-numbers, indexed by attr
        skipped = defaultdict(list)
        for run_num, r in enumerate(reductions):
            for attr in attrs:
                if attr in r:
                    accum[attr].append(r[attr])
                else:
                    accum[attr].append(place_holder) # fill in later
                    skipped[attr].append(run_num)
        if len(skipped) > 0:
            skipped_runs = list(reduce(lambda runs, runs_new: runs.union(set(runs_new)), six.itervalues(skipped), set()))
            print('Warning: Some attributes (%d of %d) were not present in some runs (%d of %d).' % (
                len(skipped), len(attrs), len(skipped_runs), len(reductions)))
            if VERBOSITY:
                print('Attrs = ' + ', '.join(skipped))
        for attr in skipped:
            # find a non-placeholder value in that attribute's values
            valid_value = None
            for result in accum[attr]:
                if result is not place_holder:
                    valid_value = result
                    break
            assert valid_value is not None, 'Attribute "%s" skipped in all runs' % attr
            # replace the place-holders with a value of the proper shape (we have naked floats and np.ndarray's)
            null_value = self.get_attr_null_value(attr, valid_value)
            #   NB: this works for all skippable values except SNR, but the "null" value may not
            #   always be 0 * v where v is a valid value.  Note in particular that
            for run_num in skipped[attr]:
                assert accum[attr][run_num] is place_holder, 'Empty-attribute fill-in logic error'
                accum[attr][run_num] = null_value
        # 2: take means, std's, quantiles of the various attributes
        #    some QOIs can have NaNs (eg, char_snr => from empty RpL bins;
        #    h_star_det_earth_frac => from stars without Earths)
        summary = {}
        #np.set_printoptions(precision=3) # for interactive debugging
        for attr in attrs:
            # shortcut for attributes that are not collapsed with a mean, std, etc.
            if attr.endswith('_list'):
                # flatten list-of-lists
                summary[attr] = [item for sublist in accum[attr] for item in sublist]
                continue
            # suppress "empty slice" warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                summary[attr + '_mean'] = np.nanmean(accum[attr], axis=0)
                # ddof=1: we're estimating the mean separately, and want Nobs=1 => no std
                summary[attr + '_std']  = np.nanstd( accum[attr], axis=0, ddof=1)
                # number of non-NaN entries in each bin of the above averages (a vector)
                N_valid = np.sum(np.isfinite(accum[attr]), axis=0)
                summary[attr + '_nEns'] = N_valid
                # sem = standard error of the _mean above
                # this is generally a vector divided by a vector; in case of 0 or 1
                # observation it will (correctly) be NaN/0 or NaN/1 => NaN.
                summary[attr + '_sem']  = summary[attr + '_std'] / np.sqrt(N_valid)
                # find q, then stop suppressing RuntimeWarnings
                # this emits a warning when a given slice is all-nan, which is not a problem
                q = np.nanpercentile(accum[attr], [25,50,75], axis=0) # skips NaN
            # scalar attributes (exoE_char_snr) will yield 1d results,
            # and 2d indexing will error out - fix it with try/except
            try:
                summary[attr + '_q25'] = q[0,:]
                summary[attr + '_q50'] = q[1,:] # median
                summary[attr + '_q75'] = q[2,:]
            except IndexError:
                if len(q) == 0:
                    q = np.zeros(3) + np.nan # happens for 0-sim runs
                summary[attr + '_q25'] = q[0]
                summary[attr + '_q50'] = q[1] # median
                summary[attr + '_q75'] = q[2]
        
        # unique detections and characterizations are handled separately
        # TODO: may want per-band chars
        if len(reductions) > 0:
            summary['dets_unique_mean'] = (1.0 *
                            sum([len(r['dets_unique' ]) for r in reductions])) / len(reductions)
            summary['chars_unique_mean'] = (1.0 *
                            sum([len(r['chars_unique']) for r in reductions])) / len(reductions)
            summary['chars_strict_mean'] = (1.0 *
                            sum([len(r['chars_strict']) for r in reductions])) / len(reductions)
        else:
            summary['dets_unique_mean' ] = 0
            summary['chars_unique_mean'] = 0
            summary['chars_strict_mean'] = 0

        # record these summaries in the object
        self.summary = summary
        self.char_bands_seen = char_bands_seen

    def dump(self, args):
        r'''Dump reduced data to files.'''

        def ensure_permissions(fn):
            r'''Ensure correct permissions on the named data file.  
            We use rw-rw-r-- = 664 (octal), to allow group-write.'''
            try:
                os.chmod(fn, 0o664)
            except OSError:
                pass # e.g., don't own the file

        # Implementation note:
        # The scheme used for dumping to each file is to define a list of QOIs, or
        # "quantities of interest", such as "exo-Earth full characterizations".
        # Then, a list of metrics-to-be-dumped is made by appending a series
        # of metrics ("_mean", "_std") to each QOI.  The CSV contains one column
        # for each QOI + metric combination.

        # basic summary statistics: mean, standard deviation, standard error of the mean
        # number of items in the cross-ensemble average
        basic_stats = ('mean', 'std', 'sem', 'nEns')

        # 0: metadata and scalars
        fn = args.outfile % ('info', 'csv')
        print('\tDumping to %s' % fn)
        # time of newest simulation (modtime of drm directory)
        simtime = '2000-01-02_03:04'
        if len(args.infile) > 0:
            drm0 = os.path.dirname(args.infile[0])
            if os.path.isdir(drm0):
                simtime = time.strftime("%Y-%m-%d_%H:%M", time.localtime(os.path.getmtime(drm0)))
        info = dict(user=os.environ['USER'],
                    runtime=time.strftime("%Y-%m-%d_%H:%M"),
                    simtime=simtime,
                    experiment=args.expt_name,
                    ensemble_size=self.Ndrm_actual,
                    detections_unique_mean=self.summary['dets_unique_mean'],
                    chars_unique_mean=self.summary['chars_unique_mean'],
                    chars_strict_mean=self.summary['chars_strict_mean'],
                    detections_earth_unique=(self.summary['exoE_det_main_mean']  +
                                             self.summary['exoE_det_alt_mean']),
                    detections_earth_all=(self.summary['exoE_xdet_main_mean'] +
                                          self.summary['exoE_xdet_alt_mean']),
                    chars_earth_unique=(self.summary['exoE_char_full_mean'] +
                                        self.summary['exoE_char_part_mean']),
                    chars_earth_strict=self.summary['exoE_char_strict_mean'],
                    # target depletion metrics -- if no chars, they will not exist
                    # these eventually land in a .json, and NaNs in json aren't OK
                    targ_dep_slope_all=self.summary.get(
                        'tdep_slope_ym1_char_full_allplan_uniq_union_mean', 0.0),
                    targ_dep_t80_all=self.summary.get(
                        'tdep_t80_char_full_allplan_uniq_union_mean', 0.0),
                        )
        with open(fn, 'w') as csvfile:
            w = csv.DictWriter(csvfile, fieldnames=sorted(info.keys()))
            w.writeheader()
            w.writerow(info)
        ensure_permissions(fn)

        # 0b: exo-Earth (scalars)
        fn = args.outfile % ('earth', 'csv')
        print('\tDumping to %s' % fn)
        # stuff to appear first
        # TODO: might want an overall count here (det_main + det_alt)?
        #       (for now, trust downstream code to add them itself)
        earth_fields = []
        # qoi = quantities-of-interest
        # compose the names of more fields to save -- originally these had no quantiles
        #   these are created by the analysis() routine and relayed to here
        #   (we modify it, so make a copy)
        earth_qoi = copy.copy(self.auto_keys.get('earth', []))
        # to this, append per-band characterization QOIs
        for band in self.char_bands_seen:
            earth_qoi.append('exoE_char_full_' + band)
            earth_qoi.append('exoE_char_part_' + band)
            earth_qoi.append('exoE_char_snr_'  + band)
        # supply quantiles for everything
        earth_fields.extend([ (x + '_' + y)
                              for x in earth_qoi
                              for y in basic_stats + ('q25', 'q50', 'q75')])
        with open(fn, 'w') as csvfile:
            w = csv.DictWriter(csvfile, fieldnames=earth_fields)
            w.writeheader()
            # dictionary mapping field -> value -- everything is a scalar here
            d = {f:self.summary[f] for f in earth_fields}
            w.writerow(d)
        ensure_permissions(fn)

        # 1: yield and resources vs. mission elapsed time
        fn = args.outfile % ('times', 'csv')
        print('\tDumping to %s' % fn)
        # make the bins appear first
        time_fields = ['h_det_time_lo', 'h_det_time_hi']
        # qoi = quantities-of-interest
        time_qoi = (
            self.auto_keys.get('delta_v', []) +
            self.auto_keys.get('resource_analysis', []))
        # compose the names of all other fields to be saved
        time_fields.extend([ (x + '_' + y)
                           for x in time_qoi
                           for y in basic_stats])
        time_len = len(self.summary[time_fields[0]])
        with open(fn, 'w') as csvfile:
            w = csv.DictWriter(csvfile, fieldnames=time_fields)
            w.writeheader()
            for i in range(time_len):
                # dictionary mapping field -> value
                d = {f:self.summary[f][i] for f in time_fields}
                w.writerow(d)
        ensure_permissions(fn)

        # 2: radius-luminosity histograms
        fn = args.outfile % ('radlum', 'csv')
        print('\tDumping to %s' % fn)
        # start making list of all fields-names to save
        #  make the bins appear first
        radlum_fields = ['h_Rp_lo', 'h_Rp_hi', 'h_L_lo', 'h_L_hi']
        # compose the names of more fields to save -- originally these had no quantiles
        #   these are created by the analysis() routine and relayed to here
        #   (we modify it, so make a copy)
        radlum_qoi = copy.copy(self.auto_keys.get('radlum', []))
        # add more names
        for band in self.char_bands_seen:
            radlum_qoi.append('h_RpL_char_full_' + band)
            radlum_qoi.append('h_RpL_char_part_' + band)
        radlum_fields.extend([ (x + '_' + y)
                              for x in radlum_qoi
                              for y in basic_stats + ('q25', 'q50', 'q75')])
        # yet more field names - these have quantiles
        radlum_qoi = []
        for band in self.char_bands_seen:
            radlum_qoi.append('h_RpL_char_snr_' + band)
        radlum_fields.extend([ (x + '_' + y)
                              for x in radlum_qoi
                              for y in basic_stats + ('q25', 'q50', 'q75')])
        # write them
        radlum_len = len(self.summary[radlum_fields[0]])
        with open(fn, 'w') as csvfile:
            w = csv.DictWriter(csvfile, fieldnames=radlum_fields)
            w.writeheader()
            for i in range(radlum_len):
                # dictionary mapping field -> value
                d = {f:self.summary[f][i] for f in radlum_fields}
                w.writerow(d)
        ensure_permissions(fn)

        # 2.5: Earth-char-count histograms
        fn = args.outfile % ('earth-char-count', 'csv')
        print('\tDumping to %s' % fn)
        # make list of all fields-names to save
        #  make the bins appear first
        earth_char_count_fields = ['h_earth_char_count_lo', 'h_earth_char_count_hi']
        # more fields to save
        earth_char_count_qoi = [
            'h_earth_char_all',
            'h_earth_xchar_all',
            'h_earth_char_strict',
            ]
        earth_char_count_fields.extend([ (x + '_' + y)
                                             for x in earth_char_count_qoi
                                             for y in basic_stats])
        # write them
        earth_char_count_len = len(self.summary[earth_char_count_fields[0]])
        with open(fn, 'w') as csvfile:
            w = csv.DictWriter(csvfile, fieldnames=earth_char_count_fields)
            w.writeheader()
            for i in range(earth_char_count_len):
                # dictionary mapping field -> value
                d = {f:self.summary[f][i] for f in earth_char_count_fields}
                w.writerow(d)
        ensure_permissions(fn)

        # 3: revisits-by-stars
        fn = args.outfile % ('visits', 'csv')
        print('\tDumping to %s' % fn)
        # compose the names of all fields to be saved
        # make the bins appear first
        visit_fields = ['h_visit_bins']
        # qoi = quantities-of-interest
        visit_qoi = self.auto_keys.get('summarize_revisits', [])
        # add the names of more fields to be saved
        visit_fields.extend([ (x + '_' + y)
                              for x in visit_qoi
                              for y in basic_stats + ('q25', 'q50', 'q75')])
        visit_len = len(self.summary[visit_fields[0]])
        with open(fn, 'w') as csvfile:
            w = csv.DictWriter(csvfile, fieldnames=visit_fields)
            w.writeheader()
            for i in range(visit_len):
                # dictionary mapping field -> value
                d = {f:self.summary[f][i] for f in visit_fields}
                w.writerow(d)
        ensure_permissions(fn)

        # 4: per-star yield, t_int, etc.
        fn = args.outfile % ('star-target', 'csv')
        print('\tDumping to %s' % fn)
        # compose the names of all fields to be saved
        # make the names and other star properties appear first
        star_target_fields = ['star_sInd', 'star_Name', 'star_L', 'star_dist', 'star_Spec']
        # qoi = quantities-of-interest
        # set up in two routines, per_star_yield()/..._promotion(), and relayed to here
        star_target_qoi = (
            self.auto_keys.get('per_star_yield', []) +
            self.auto_keys.get('per_star_promotion', []))
        # add the names of more fields to be saved
        # (now using std in javascript plot widget)
        star_target_fields.extend([ (x + '_' + y)
                              for x in star_target_qoi
                              for y in basic_stats])
        star_target_len = len(self.summary[star_target_fields[0]])
        with open(fn, 'w') as csvfile:
            w = csv.DictWriter(csvfile, fieldnames=star_target_fields)
            w.writeheader()
            for i in range(star_target_len):
                # dictionary mapping field -> value
                d = {f:self.summary[f][i] for f in star_target_fields}
                w.writerow(d)
        ensure_permissions(fn)

        # 5: event durations
        fn = args.outfile % ('events', 'csv')
        print('\tDumping to %s' % fn)
        # compose the names of all fields to be saved
        # make the bins appear first
        event_fields = ['h_event_b0_duration_lo', 'h_event_b0_duration_hi',
                        'h_event_b1_duration_lo', 'h_event_b1_duration_hi',
                        'h_event_b2_duration_lo', 'h_event_b2_duration_hi']
        # qoi = quantities-of-interest
        # set up in event_analysis() and relayed to here
        event_qoi = self.auto_keys.get('event_analysis', [])
        # add the names of more fields to be saved
        # add quantiles for durations, which may have large values
        event_fields.extend([ (x + '_' + y)
                              for x in event_qoi
                              for y in basic_stats + ('q25', 'q50', 'q75')])
        event_len = len(self.summary[event_fields[0]])
        with open(fn, 'w') as csvfile:
            w = csv.DictWriter(csvfile, fieldnames=event_fields)
            w.writeheader()
            for i in range(event_len):
                # dictionary mapping field -> value
                d = {f:self.summary[f][i] for f in event_fields}
                w.writerow(d)
        ensure_permissions(fn)

        # 5: event counts (as histograms)
        fn = args.outfile % ('event-counts', 'csv')
        print('\tDumping to %s' % fn)
        # compose the names of all fields to be saved
        # make the bins appear first
        event_count_fields = ['h_event_count_lo', 'h_event_count_hi']
        # set up in count_events() and relayed to here
        event_count_qoi = self.auto_keys.get('event_count', [])
        # add the names of more fields to be saved
        event_count_fields.extend([ (x + '_' + y)
                                for x in event_count_qoi
                                for y in basic_stats])
        event_count_len = len(self.summary[event_count_fields[0]])
        with open(fn, 'w') as csvfile:
            w = csv.DictWriter(csvfile, fieldnames=event_count_fields)
            w.writeheader()
            for i in range(event_count_len):
                # dictionary mapping field -> value
                d = {f:self.summary[f][i] for f in event_count_fields}
                w.writerow(d)
        ensure_permissions(fn)

        # 6A: planet promotion counts (each QOI is a function of time)
        fn = args.outfile % ('promote', 'csv')
        print('\tDumping to %s' % fn)
        # compose the names of all fields to be saved
        # make the bins appear first
        promo_fields = ['h_promo_time_lo', 'h_promo_time_hi']
        # qoi = quantities-of-interest
        #   these are created by the promotion_analysis() routine and relayed to here
        promo_qoi = self.auto_keys.get('promo_count', [])
        # add the names of more fields to be saved
        promo_fields.extend([ (x + '_' + y)
                                for x in promo_qoi
                                for y in basic_stats + ('q25', 'q50', 'q75')])
        promo_len = len(self.summary[promo_fields[0]])
        # allow to skip if the promotion analysis could not be done.
        if promo_fields[-1] in self.summary:
            with open(fn, 'w') as csvfile:
                w = csv.DictWriter(csvfile, fieldnames=promo_fields)
                w.writeheader()
                for i in range(promo_len):
                    # dictionary mapping field -> value
                    d = {f:self.summary[f][i] for f in promo_fields}
                    w.writerow(d)
            ensure_permissions(fn)
        else:
            print('\tDump to %s: skipped.' % fn)
            
        # 6B: planet promotion planet-count histograms
        # (each QOI is a histogram of counts, at a handful of pre-chosen time slices)
        fn = args.outfile % ('promote-hist', 'csv')
        print('\tDumping to %s' % fn)
        # compose the names of all fields to be saved
        # make the bins appear first
        phist_fields = ['h_phist_count_lo', 'h_phist_count_hi']
        # qoi = quantities-of-interest
        #   these are created by the promotion_analysis() routine and relayed to here
        phist_qoi = self.auto_keys.get('promo_phist', [])
        # add the names of more fields to be saved
        phist_fields.extend([ (x + '_' + y)
                                for x in phist_qoi
                                for y in ('mean', )])
        phist_len = len(self.summary[phist_fields[0]])
        # allow to skip if the promotion analysis could not be done.
        if phist_fields[-1] in self.summary:
            with open(fn, 'w') as csvfile:
                w = csv.DictWriter(csvfile, fieldnames=phist_fields)
                w.writeheader()
                for i in range(phist_len):
                    # dictionary mapping field -> value
                    d = {f:self.summary[f][i] for f in phist_fields}
                    w.writerow(d)
            ensure_permissions(fn)
        else:
            print('\tDump to %s: skipped.' % fn)
            
        # 7: funnel analysis (promotion, deep-dive)
        fn = args.outfile % ('funnel', 'csv')
        print('\tDumping to %s' % fn)
        # compose the names of all fields to be saved
        # initial list
        funnel_fields = []
        # qoi = quantities-of-interest
        #   these are created by the funnel_analysis() routine and relayed to here
        funnel_qoi = self.auto_keys.get('funnel', [])
        # add the names of more fields to be saved
        funnel_fields.extend([ (x + '_' + y)
                                for x in funnel_qoi
                                for y in basic_stats])
        # NB: funnel_* quantities are scalars!
        # funnel_len = len(self.summary[funnel_fields[0]])
        with open(fn, 'w') as csvfile:
            w = csv.DictWriter(csvfile, fieldnames=funnel_fields)
            w.writeheader()
            # dictionary mapping field -> value
            d = {f:self.summary[f] for f in funnel_fields}
            w.writerow(d)
        ensure_permissions(fn)

        # 7.2: det_funnel analysis (allstars, promotion)
        fn = args.outfile % ('det-funnel', 'csv')
        print('\tDumping to %s' % fn)
        # compose the names of all fields to be saved
        # initial list
        detfunnel_fields = []
        # qoi = quantities-of-interest
        #   these are created by the det_funnel_analysis() routine and relayed to here
        detfunnel_qoi = self.auto_keys.get('detfunnel', [])
        # add the names of more fields to be saved
        detfunnel_fields.extend([ (x + '_' + y)
                                for x in detfunnel_qoi
                                for y in basic_stats])
        # NB: detfunnel_* quantities are scalars!
        # funnel_len = len(self.summary[funnel_fields[0]])
        with open(fn, 'w') as csvfile:
            w = csv.DictWriter(csvfile, fieldnames=detfunnel_fields)
            w.writeheader()
            # dictionary mapping field -> value
            d = {f:self.summary[f] for f in detfunnel_fields}
            w.writerow(d)
        ensure_permissions(fn)

        # 8: earth-char-attempts analysis
        fn = args.outfile % ('earth-char-list', 'csv')
        print('\tDumping to %s' % fn)
        # named fields for earth-char list (it's just one field in self.summary, hence [0])
        earth_char_qoi = self.auto_keys.get('earth_char', [])[0]
        # data for earth-char list (it's just one field in self.summary)
        earth_char_data = self.summary[earth_char_qoi]
        # compose the names of all fields to be saved
        # initial list
        earth_char_fields = list(earth_char_data[0].keys()) if earth_char_data else []
        # NB: earth_char_* quantities are scalars!
        with open(fn, 'w') as csvfile:
            w = csv.DictWriter(csvfile, fieldnames=earth_char_fields)
            w.writeheader()
            # dictionary mapping field -> value
            for row in earth_char_data:
                # dictionary mapping field -> value
                w.writerow(row)
        ensure_permissions(fn)

        # 9: yield-vs-time analysis
        fn = args.outfile % ('yield-time', 'csv')
        print('\tDumping to %s' % fn)
        # compose the names of all fields to be saved
        # initial list, will appear first in output
        yield_time_fields = ['h_det_time_lo', 'h_det_time_hi']
        # qoi = quantities-of-interest
        #   these are created by the yield_time_analysis() routine and relayed to here
        yield_time_qoi = self.auto_keys.get('yield_time', [])
        # add the names of more fields to be saved
        yield_time_fields.extend([ (x + '_' + y)
                                for x in yield_time_qoi
                                for y in basic_stats])
        yield_time_len = len(self.summary[yield_time_fields[0]])
        with open(fn, 'w') as csvfile:
            w = csv.DictWriter(csvfile, fieldnames=yield_time_fields)
            w.writeheader()
            for i in range(time_len):
                # dictionary mapping field -> value
                d = {f:self.summary[f][i] for f in yield_time_fields}
                w.writerow(d)
        ensure_permissions(fn)

        # 10: visits-vs-time analysis
        fn = args.outfile % ('visit-time', 'csv')
        print('\tDumping to %s' % fn)
        # compose the names of all fields to be saved
        # initial list, will appear first in output
        visit_time_fields = ['h_det_time_lo', 'h_det_time_hi']
        # qoi = quantities-of-interest
        #   these are created by visit_time_analysis() and relayed to here
        visit_time_qoi = self.auto_keys.get('visit_time', [])
        # add the names of more fields to be saved
        visit_time_fields.extend([ (x + '_' + y)
                                for x in visit_time_qoi
                                for y in basic_stats])
        visit_time_len = len(self.summary[visit_time_fields[0]])
        with open(fn, 'w') as csvfile:
            w = csv.DictWriter(csvfile, fieldnames=visit_time_fields)
            w.writeheader()
            for i in range(time_len):
                # dictionary mapping field -> value
                d = {f:self.summary[f][i] for f in visit_time_fields}
                w.writerow(d)
        ensure_permissions(fn)

        # 11: target-depletion analysis
        fn = args.outfile % ('target-depletion', 'csv')
        print('\tDumping to %s' % fn)
        # compose the names of all fields to be saved
        # initial list
        target_depletion_fields = []
        # qoi = quantities-of-interest
        #   these are created by the ..._analysis() routine and relayed to here
        target_depletion_qoi = self.auto_keys.get('target_depletion', [])
        # add the names of more fields to be saved
        target_depletion_fields.extend([ (x + '_' + y)
                                    for x in target_depletion_qoi
                                    for y in basic_stats])
        # NB: target_depletion_* quantities are scalars!
        with open(fn, 'w') as csvfile:
            w = csv.DictWriter(csvfile, fieldnames=target_depletion_fields)
            w.writeheader()
            # dictionary mapping field -> value
            d = {f:self.summary[f] for f in target_depletion_fields}
            w.writerow(d)
        ensure_permissions(fn)

        # copy script to JSON file
        fn = args.outfile % ('script', 'json')
        print('\tCopying script to %s' % fn)
        shutil.copyfile(args.script_name, fn)
        ensure_permissions(fn)
                
        # copy an outspec to specific JSON file
        fn = args.outfile % ('outspec', 'json')
        print('\tCopying outspec to %s' % fn)
        shutil.copyfile(args.specs_name, fn)
        ensure_permissions(fn)

        # dump more info to JSON file
        fn = args.outfile % ('pool', 'json')
        pool = self.summary
        with open(fn, 'w') as outfile:
            json.dump(pool, outfile, sort_keys=True, indent=4, ensure_ascii=False,
                        separators=(',', ': '), default=array_encoder)
        ensure_permissions(fn)

                
def obtain_outspec(args, announce=True):
    r'''Identify a good "specs" (script) file using Sandbox structure

    If the script takes Exosims defaults, it will not contain explicit
    values for some parameters (like settlingTime). But the outspec will.
    This function tries to find the outspec within the logfiles, in
    preference to the script, so that the defaults will be explicit.'''
    # try to get some outspec nearby the .pkl, else, the script itself

    # parent of parent:
    # sims/Yokohama_Extended.fam/Trials.fam/s_dbug_YX_NIR_D8.0/drm/397016230.pkl
    #   -> 
    # sims/Yokohama_Extended.fam/Trials.fam/s_dbug_YX_NIR_D8.0
    sim_root = os.path.dirname(os.path.dirname(args.infile[0]))
    if os.path.isdir(spec_dir := os.path.join(sim_root, 'run')):
        # 2022'ish logfiles: .../run/outspec_SEED.json
        path_pat = os.path.join(spec_dir, 'outspec_*.json')
    elif os.path.isdir(spec_dir := os.path.join(sim_root, 'log', 'outspec')):
        # 2023+ logfiles: .../log/outspec/SEED.json
        path_pat = os.path.join(spec_dir, '*.json')
    else:
        path_pat = ''
    # settle for anything
    paths = glob.glob(path_pat)
    # try to do better than the original script
    rv = paths[0] if paths else args.script_name
    if paths and announce:
        print(f"{args.progname}: Params from outspec: {rv}")
    return rv


def load_exosims_sim(args):
    r'''Extract fields from an Exosims sim corresponding to the DRMs being reduced.

    We return a "god object" with only the fields we need present.
    For now: ohTime and settlingTime.
    None is returned if the script could not be loaded. '''

    # Developer note: the "specs" may be the original (human-compiled) input script,
    # or the outspec as written by Exosims. The code here must tolerate either.

    # identify and load a good "specs"
    args.specs_name = obtain_outspec(args, announce=True)
    try:
        with open(args.specs_name, 'r') as fp:
            specs = json.load(fp)
    except IOError:
        print('%s: Could not open implied script/specs file <%s>.' % (args.progname, args.specs_name))
        print('%s: Promotion analysis will be skipped.', (args.progname, ))
        return None

    # formulate the simulation-summary object -- only the fields we need
    rv = {}

    # 1: deep-dive targets => 'top_HIPs'
    if 'occHIPs' not in specs:
        top_HIPs = []
    else:
        # obtain occHIPs, the list of Hipparcos numbers
        if isinstance(specs['occHIPs'], list):
            # it is a literal list
            occHIPs = specs['occHIPs']
        else:
            # it is a filename
            occHIPs_path = specs['occHIPs']
            # occHIPs_path = os.path.join(EXOSIMS.__path__[0], 'Scripts', occHIPs_path)
            try:
                with open(occHIPs_path, 'r') as ofile:
                    HIPs_text = ofile.read()
            except IOError:
                print("%s: Could not open occHIPs file `%s'." % (args.progname, occHIPs_path))
                print("%s: Continuing with null occHIPs list." % (args.progname, ))
                HIPs_text = ''
            if ',' in HIPs_text:
                occHIPs = HIPs_text.split(',')
            else:
                occHIPs = HIPs_text.split('\n')
        # strip surrounding spaces
        occHIPs = [hip.strip() for hip in occHIPs]
        # restrict to topstars
        topstars = specs.get('topstars', 0)
        top_HIPs = occHIPs[:topstars]
    rv['top_HIPs'] = top_HIPs

    # 2: obtain overhead time
    # note: we do not strictly need to instantiate the object to get this
    # note: the module instantiation set up the ohTime as an astropy Quantity
    # note: this is actually a function of the SuppressionSystem
    ohTime = 1.0 # [days] --- this is the Exosims default
    for syst in specs.get('starlightSuppressionSystems', []):
        ohTime = max(ohTime, strip_units(syst.get('ohTime', 0.0)))
    rv['ohTime'] = ohTime
    # obtain settling time
    # settling time is in Observatory
    #  -- 1.0 days is the Exosims default
    settlingTime = specs.get('settlingTime', 1.0) # [days]
    rv['settlingTime'] = settlingTime

    # 3: obtain mission lifetime
    rv['missionLife'] = specs.get('missionLife', 5.0) # [years]

    return rv

def main(args):
    print('%s: Loading %d file patterns.' % (args.progname, len(args.infile)))
    args.sim_info = load_exosims_sim(args)
    UPDATE_GLOBALS(args.sim_info)
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
    parser.add_argument('drm', metavar='DRM', nargs='*', default=[],
                            help='drm file(s)')
    parser.add_argument('-s', '--script', type=str, default='',
                            help='Exosims JSON script used for these DRMs.')
    parser.add_argument('-O', '--outfile', type=str, default='',
                            help='Output file template.')
    parser.add_argument('-v', help='verbosity, repeat to escalate', action='count', dest='verbose',
                            default=0)
    parser.add_argument('-D', help='debugging output type (string)', type=str, default='', dest='debug')
    parser.add_argument('-j', '--jobs', help=' #jobs, default = %(default)d',
                      type=int, dest='jobs', default=max(1,int(N_CPU*0.35)))
    args = parser.parse_args()
    
    # measure elapsed time
    time_start = time.time()

    # set umask in hopes that files will be group-writable
    os.umask(0o002)

    VERBOSITY = args.verbose
    DEBUG = args.debug

    # program name, for convenience
    args.progname = os.path.basename(sys.argv[0])
    
    # preserve input args separately
    args.infile = args.drm[:]

    # check edge cases
    if len(args.infile) == 1 and '*' in os.path.basename(args.infile[0]):
        # the wildcard appears here when invoked with arg .../drm/*.pkl, and there are no .pkl's
        # this fix is a bit hacky (there "could" be a drm with a star in its name),
        # but it helps do the right thing with automatic invocation by make
        print('%s: Wildcard appears in input drm, assuming empty.' % args.progname)
        print('%s: Subsequent reduction error is likely.' % args.progname)
        args.infile = []
    if len(args.infile) == 0:
        print('%s: No input DRMs. Exiting.' % args.progname)
        sys.exit(0)

    # get the experiment name from the directory
    #   this is brittle, but have to do something
    try:
        # ' ' + == hack to distinguish auto-generated from manual expt_name
        args.expt_name = ' ' + os.path.dirname(args.drm[0]).split('/')[-2]
    except IndexError:
        args.expt_name = ' EXOSIMS Run'
    # allow to over-ride with a single name in a special file
    try:
        ens_fn = re.sub(r'/drm/.*', '/EnsembleName.txt', args.drm[0])
        if os.path.exists(ens_fn):
            args.expt_name = open(ens_fn).readline().strip()
            print('%s: New ensemble name is "%s".' % (args.progname, args.expt_name))
    except:
        pass

    # get the script name from the directory - this relies heavily on Sandbox convention
    # TODO: allow to specify manually?
    if not args.script:
        args.script_name = os.path.dirname(args.drm[0]).replace('sims/', 'Scripts/').replace('/drm', '.json')
        print('%s: Inferring script file: %s' % (args.progname, args.script_name))
    else:
        args.script_name = args.script
        print('%s: Explicitly given script file: %s' % (args.progname, args.script_name))

    # best practice is to explicitly give outfile
    if not args.outfile:
        args.outfile = ('sims/%s/reduce' % args.expt_name) + '-%s.%s'
    if args.outfile.count('%s') != 2:
        print('%s: Need two instances of %%s in output file template' % args.progname)
        sys.exit(1)
    # ensure enclosing dir exists
    directory = os.path.dirname(args.outfile % ('dummy', 'txt'))
    if not os.path.exists(directory):
        os.makedirs(directory)

    # if this file is present, do not continue
    skip_fn = args.outfile % ('skip', 'txt')
    if os.path.isfile(skip_fn):
        print('%s: Skipping reduction in %s.' % (args.progname, os.path.dirname(args.outfile)))
        sys.exit(0)

    infile_print = (args.infile[0] if len(args.infile) > 0 else '(none)') + (' ...' if len(args.infile) > 1 else '')
    print('%s: Reducing %s to %s' % (args.progname, infile_print, args.outfile))

    main(args)
    # measure elapsed time
    time_delta = time.time() - time_start
    print(f'{args.progname}: Done. Elapsed time: {time_delta:.3}s.')

    sys.exit(0)
