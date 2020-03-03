#!/usr/bin/env python
"""plot-timeline: plot timeline of all observations within a single DRM

Usage:

  plot-timeline.py [-d] [-s] [-t CSV] [-o OUTPATH] SCRIPT DRM

where:

  -d means to omit the day-in-the-life plot (can be useful)
  -s means to omit the synoptic all-mission plot
  -o OUTPATH gives an output path for results, containing two %s slots  
  -t CSV names a CSV file for General Astrophysics proxy observation-durations
     (without -t, a sane builtin default is used).

Most helpful Sandbox usage:

  plot-timeline.py -d Scripts/FOO.json sims/FOO/drm/SEED.pkl 

where FOO.json is a script, and SEED.pkl is a DRM.  Output will be placed in the working directory 
unless -o path/to/output/%s.%s or the like is given.

More-complex actual usage line:

  util/plot-timeline.py -d Scripts/HabEx_4m_TSDDold_DD_TF17_maxchar3_a0.78b0.2c0.07d0.11e.58f.05__20190402.json sims/HabEx_4m_TSDDold_DD_TF17_maxchar3_a0.78b0.2c0.07d0.11e.58f.05__20190402/drm/992454934.pkl

Michael Turmon, JPL, 04/2019 -- created based on an idea by Dean Keithly

"""

from __future__ import print_function
from six.moves import range
from six.moves import zip
import six.moves.cPickle as pickle
import os
import sys
import csv
import argparse
import json
import random
import time
import numpy as np

# currently this pylab import is needed to allow the script
# to terminate cleanly
#from pylab import *
#import matplotlib
#matplotlib.use('Agg')
import matplotlib as mpl; mpl.use('Agg') # not interactive: don't use X backend
import matplotlib.pyplot as plt

# unpickling python2/numpy pickles within python3 requires this
PICKLE_ARGS = {} if sys.version_info.major < 3 else {'encoding': 'latin1'}


########################################
###
###  Static Data
###
########################################

# This is a histogram of observation-durations that is sampled to fill in unknown times for General Astrophysics.
# Use the "-t" option to specify a CSV file of your choosing, to replace these times, if you wish.
# The times here were compiled from a list of Spitzer observation durations provided by Rashied Amini of JPL, covering
# 2009-2018, using the script "spitz-parse.py", which outputs a CSV file containing columns:
#   (left bin edge, in minutes) (probability mass within bin) (CDF at right edge of bin)
# This CSV file (~200 lines) was condensed into the values below by this command:
#   awk -F, '{printf("(%4s,%11s, %7s),\n", $1, $2, $3)}' spitzer-obs-dur-hist.csv | paste -d " " - - - - | sed 's/^/    /'
# Column 3 is the runing sum of column 2, and we only use columns 1 and 2.
# Doing this allows this script to be self-contained.

DEFAULT_OBS_DURATIONS = [
    (   4,    0.56698, 0.56698), (  19,    0.16423, 0.73121), (  34,   0.049178, 0.78039), (  49,   0.040588, 0.82098),
    (  64,   0.032185, 0.85316), (  79,   0.018549, 0.87171), (  94,   0.034626, 0.90634), ( 109,   0.022414, 0.92875),
    ( 124,  0.0057415, 0.93449), ( 139,  0.0046153, 0.93911), ( 154,  0.0036547, 0.94276), ( 169,  0.0042951, 0.94706),
    ( 184,  0.0026278, 0.94968), ( 199,  0.0031136,  0.9528), ( 214, 0.00085018, 0.95365), ( 229, 0.00093851, 0.95459),
    ( 244,  0.0043613, 0.95895), ( 259,   0.004284, 0.96323), ( 274,  0.0049134, 0.96815), ( 289,  0.0018991, 0.97004),
    ( 304,  0.0034338, 0.97348), ( 319,  0.0028707, 0.97635), ( 334,  0.0021199, 0.97847), ( 349, 0.00043061,  0.9789),
    ( 364, 0.00076185, 0.97966), ( 379,  0.0027272, 0.98239), ( 394, 0.00059623, 0.98299), ( 409, 0.00034228, 0.98333),
    ( 424, 0.00086122, 0.98419), ( 439,  0.0005079,  0.9847), ( 454, 0.00061831, 0.98532), ( 469, 0.00070664, 0.98602),
    ( 484, 0.00086122, 0.98688), ( 499, 0.00041957,  0.9873), ( 514, 0.00076185, 0.98806), ( 529, 0.00071768, 0.98878),
    ( 544, 0.00086122, 0.98964), ( 559, 0.00023187, 0.98988), ( 574, 0.00025395, 0.99013), ( 589, 0.00054102, 0.99067),
    ( 604, 0.00083914, 0.99151), ( 619,  0.0003202, 0.99183), ( 634, 0.00015458, 0.99198), ( 649, 0.00019874, 0.99218),
    ( 664, 0.00015458, 0.99234), ( 679, 0.00025395, 0.99259), ( 694, 0.00016562, 0.99276), ( 709, 0.00017666, 0.99293),
    ( 724, 0.00038645, 0.99332), ( 739, 0.00020978, 0.99353), ( 754, 0.00011041, 0.99364), ( 769, 0.00016562, 0.99381),
    ( 784, 6.6248e-05, 0.99387), ( 799,  0.0001325,   0.994), ( 814, 3.3124e-05, 0.99404), ( 829, 0.00061831, 0.99466),
    ( 844, 0.00017666, 0.99483), ( 859, 6.6248e-05,  0.9949), ( 874, 9.9372e-05,   0.995), ( 889, 0.00023187, 0.99523),
    ( 904, 0.00012145, 0.99535), ( 919, 3.3124e-05, 0.99538), ( 934, 0.00015458, 0.99554), ( 949, 6.6248e-05, 0.99561),
    ( 964, 9.9372e-05,  0.9957), ( 979, 0.00020978, 0.99591), ( 994, 0.00015458, 0.99607), (1009, 6.6248e-05, 0.99614),
    (1024,  0.0001325, 0.99627), (1039, 3.3124e-05,  0.9963), (1054, 9.9372e-05,  0.9964), (1069, 7.7289e-05, 0.99648),
    (1084,  8.833e-05, 0.99657), (1099, 6.6248e-05, 0.99663), (1114, 4.4165e-05, 0.99668), (1129, 1.1041e-05, 0.99669),
    (1144, 3.3124e-05, 0.99672), (1159, 4.4165e-05, 0.99676), (1174, 4.4165e-05, 0.99681), (1189,  0.0001877,   0.997),
    (1204, 4.4165e-05, 0.99704), (1219, 4.4165e-05, 0.99709), (1234, 0.00033124, 0.99742), (1249,  0.0001325, 0.99755),
    (1264,  8.833e-05, 0.99764), (1279, 1.1041e-05, 0.99765), (1294, 0.00027603, 0.99792), (1309, 5.5207e-05, 0.99798),
    (1324, 2.2083e-05,   0.998), (1339, 2.2083e-05, 0.99802), (1354, 5.5207e-05, 0.99808), (1369, 1.1041e-05, 0.99809),
    (1384,          0, 0.99809), (1399, 2.2083e-05, 0.99811), (1414, 2.2083e-05, 0.99813), (1429, 6.6248e-05,  0.9982),
    (1444, 0.00019874,  0.9984), (1459, 2.2083e-05, 0.99842), (1474, 2.2083e-05, 0.99844), (1489,  0.0003754, 0.99882),
    (1504,          0, 0.99882), (1519,          0, 0.99882), (1534,          0, 0.99882), (1549, 1.1041e-05, 0.99883),
    (1564, 5.5207e-05, 0.99888), (1579,          0, 0.99888), (1594,          0, 0.99888), (1609, 2.2083e-05, 0.99891),
    (1624, 1.1041e-05, 0.99892), (1639,          0, 0.99892), (1654,          0, 0.99892), (1669,          0, 0.99892),
    (1684, 1.1041e-05, 0.99893), (1699, 2.2083e-05, 0.99895), (1714,          0, 0.99895), (1729,          0, 0.99895),
    (1744, 1.1041e-05, 0.99896), (1759,          0, 0.99896), (1774,          0, 0.99896), (1789,          0, 0.99896),
    (1804,          0, 0.99896), (1819, 1.1041e-05, 0.99897), (1834, 2.2083e-05,   0.999), (1849,          0,   0.999),
    (1864,          0,   0.999), (1879,          0,   0.999), (1894,          0,   0.999), (1909,          0,   0.999),
    (1924, 4.4165e-05, 0.99904), (1939,          0, 0.99904), (1954,          0, 0.99904), (1969,          0, 0.99904),
    (1984,          0, 0.99904), (1999,          0, 0.99904), (2014, 2.2083e-05, 0.99906), (2029,          0, 0.99906),
    (2044, 5.5207e-05, 0.99912), (2059, 4.4165e-05, 0.99916), (2074,          0, 0.99916), (2089,          0, 0.99916),
    (2104, 2.2083e-05, 0.99918), (2119,          0, 0.99918), (2134,          0, 0.99918), (2149,          0, 0.99918),
    (2164, 1.1041e-05, 0.99919), (2179, 1.1041e-05, 0.99921), (2194, 2.2083e-05, 0.99923), (2209,          0, 0.99923),
    (2224,          0, 0.99923), (2239, 4.4165e-05, 0.99927), (2254,          0, 0.99927), (2269, 1.1041e-05, 0.99928),
    (2284,          0, 0.99928), (2299,          0, 0.99928), (2314,          0, 0.99928), (2329, 1.1041e-05, 0.99929),
    (2344, 4.4165e-05, 0.99934), (2359, 1.1041e-05, 0.99935), (2374,          0, 0.99935), (2389, 2.2083e-05, 0.99937),
    (2404, 3.3124e-05,  0.9994), (2419, 2.2083e-05, 0.99943), (2434, 2.2083e-05, 0.99945), (2449,          0, 0.99945),
    (2464, 1.1041e-05, 0.99946), (2479,          0, 0.99946), (2494,          0, 0.99946), (2509, 2.2083e-05, 0.99948),
    (2524,          0, 0.99948), (2539, 1.1041e-05, 0.99949), (2554, 1.1041e-05,  0.9995), (2569,          0,  0.9995),
    (2584,          0,  0.9995), (2599,          0,  0.9995), (2614,          0,  0.9995), (2629,          0,  0.9995),
    (2644,          0,  0.9995), (2659,          0,  0.9995), (2674,          0,  0.9995), (2689, 1.1041e-05, 0.99951),
    (2704,          0, 0.99951), (2719, 2.2083e-05, 0.99954), (2734,          0, 0.99954), (2749, 1.1041e-05, 0.99955),
    (2764,          0, 0.99955), (2779,          0, 0.99955), (2794,          0, 0.99955), (2809,          0, 0.99955),
    (2824, 1.1041e-05, 0.99956), (2839,          0, 0.99956), (2854,          0, 0.99956), (2869, 0.00044165,       1),
    ]


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

def load_dt_generator_info(obstime_csv):
    r'''Load info to initialize a generator of random observation times from a CSV file.'''
    values, probs = [], []
    if obstime_csv:
        # external file
        with open(obstime_csv, 'r') as f:
            reader = csv.reader(f)
            for row in reader:
                values.append(float(row[0]))
                probs .append(float(row[1]))
    else:
        # builtin values
        for row in DEFAULT_OBS_DURATIONS:
            values.append(float(row[0]))
            probs .append(float(row[1]))
    probs_a = np.array(probs)
    # ensure unit norm on the probabilities
    # convert minutes in CSV values to days units
    return np.array(values)/(24.0*60.0), probs_a/np.sum(probs_a)

def draw_obstime(values, probs):
    r'''Return a random observation time.'''
    # random offset within the binwidth of values - assume uniform bins
    dt = values[1] - values[0]
    dither = np.random.uniform() * dt
    # choose a bin from values, and add the within-bin offset
    return np.random.choice(a=values, p=probs) + dither

########################################
###
###  Simulation Container Class
###
########################################

class SimulationRun(object):
    def __init__(self, drmfile, script):
        """ loads pkl and script files
        Args:
            drmfile (string) - path to drm pickle file to load
            script (string) - path to Exosims script .json file
        Return:
            initialized SimulationRun object

        Instantiates:
            DRM (list) - a list of observations
            outspec (dict) - a dict containing input instructions
        """
        try:
            with open(drmfile, 'rb') as f:
                DRM = pickle.load(f, **PICKLE_ARGS)
        except:
            sys.stderr.write('Failed to open DRM file "%s"\n' % drmfile)
            raise
        try:
            with open(script, 'rb') as g:
                outspec = json.load(g)
        except:
            sys.stderr.write('Failed to open script file "%s"\n' % script)
            raise
        # save this in the object state
        self.drm = DRM
        self.outspec = outspec
        self.name = os.path.splitext(os.path.basename(drmfile))[0] 
        # expedient way to suppress the numerical seed for presentation plots
        self.show_seed = True
        try:
            if os.path.exists(os.path.join(os.path.dirname(drmfile), '..', 'EnsembleName.txt')):
                self.show_seed = False
        except:
            pass
        #
        # set up auxiliary quantities
        #
        mode_det = [mode for mode in outspec['observingModes'] if 'detection' in list(mode.keys())]
        if len(mode_det) > 1:
            t_mult = mode_det[0].get('timeMultiplier', 1.0) # exosims default = 1.0
        else:
            # if no detection modes, give a warning, but continue
            t_mult = 1.0
            sys.stderr.write('No detection modes found, using timeMultiplier = %f\n' % t_mult)
        self.timeMultiplier = t_mult
        # FIXME: there is also a char_margin parameter
        # FIXME: for ohTime, use what is in the starlight suppression system
        #   (outspec['starlightSuppressionSystems'][0]['ohTime']  or so)
        self.ohTime = self.outspec.get('ohTime', 0.2) # [days] -- overhead time
        self.settlingTime = self.outspec['settlingTime'] # [days] -- settling time
        self.missionLife = self.outspec['missionLife'] * 365.25 # [days]
        self.charMargin = self.outspec.get('charMargin', 0.15) # [dimensionless] - Exosims default

    def get_obs_times(self):
        r'''Extract observation start and end times, place in object.

        Note on ohTime, settlingTime, and friends:
        - ohTime and settlingTime are waited-out at the *start* of the observing window.
        The sequence for a series of detections is:
           ohTime, settlingTime, det_time, ohTime, settlingTime, det_time, ...
        In particular, ohTime and settlingTime are *not* part of det_time
        - timeMultiplier is an integration time multiplier, equal to the number of discrete 
          integrations needed to cover the full field of view, or the full wavelength band.
        '''

        ohTime = self.ohTime
        settlingTime = self.settlingTime
        timeMultiplier = self.timeMultiplier
        # initial values
        det_t0,  det_dt  = [], []
        char_t0, char_dt = [], []
        slew_t0, slew_dt = [], []
        for obs in self.drm:
            arrival_time = strip_units(obs['arrival_time'])
            if ('det_info' in obs) or ('det_time' in obs):
                # a detection
                obs_time = strip_units(obs['det_time'])*timeMultiplier + ohTime + settlingTime
                det_t0.append(arrival_time)
                det_dt.append(obs_time)
            if has_char_info(obs):
                # a characterization -- formerly, just an "else" -- assuming that all obs
                #   were detections or characterizations
                # a characterization
                # Note: charMargin is already included in char_time
                # TODO: check if timeMultiplier is included in char_time (not critical, we don't use it)
                # TODO: it can be that char_time = 0, perhaps not include such obs at all?
                char_t0.append(arrival_time)
                char_dt.append(strip_units(get_char_time(obs)))
            # put the slew info in -- will be missing from coro-only DRMs, obs.get() covers this
            slew_time = strip_units(obs.get('slew_time', 0.0)) # [days]
            if slew_time > 0:
                # e.g., first slew is length 0
                slew_t0.append(arrival_time-slew_time)
                slew_dt.append(slew_time)
        
        self.det_t0  = np.array(det_t0)
        self.det_dt  = np.array(det_dt)
        self.char_t0 = np.array(char_t0)
        self.char_dt = np.array(char_dt)
        self.slew_t0 = np.array(slew_t0)
        self.slew_dt = np.array(slew_dt)
        # extra field for synoptic plots: round up char dt to make each one visible in plot
        self.char_dt_pad = np.maximum(self.char_dt, 0.25)
        # full mission timeline
        self.all_t0  = np.array([0.0])
        self.all_dt  = np.array([self.missionLife]) # is already in days

    def find_ga_times(self):
        r'''Find General Astrophysics start and end times, place in object.
        The GA start and end times are simply times when we are neither doing 
        detection or characterization.'''
        # the maximum det or char time window [days]
        t_final = np.maximum((self.char_t0[-1]+self.char_dt[-1]) if len(self.char_t0) > 0 else 0.0,
                             (self.det_t0 [-1]+self.det_dt [-1]) if len(self.det_t0 ) > 0 else 0.0) + 1
        # list of sampled times from 0..t_final, in fixed increments
        t_range = np.arange(0.0, t_final, 1.0/(24.0*60.0))
        # boolean for where GA observations may be taken, aligned with t_range
        ga_ok = np.full(len(t_range), True)
        # set GA *not* ok during all dets and chars
        for i in range(len(self.char_t0)):
            t0 = self.char_t0[i]
            t1 = self.char_t0[i] + self.char_dt[i]
            ga_ok[np.logical_and(t_range >= t0, t_range <= t1)] = False
        for i in range(len(self.det_t0)):
            t0 = self.det_t0[i]
            t1 = self.det_t0[i] + self.det_dt[i]
            ga_ok[np.logical_and(t_range >= t0, t_range <= t1)] = False
        # ga_ok goes [False True...True False] and np.diff picks out edges
        # ensure we start and end at False, so that the "diff" method works
        ga_ok[0] = False
        ga_ok[-1] = False
        breaks = np.diff(ga_ok.astype(np.int16))
        ga_t0 = t_range[np.nonzero(breaks ==  1)]
        ga_dt = t_range[np.nonzero(breaks == -1)] - ga_t0
        if False:
            # use slew time windows instead
            self.ga_t0 = self.slew_t0
            self.ga_dt = self.slew_dt
        else:
            # preserve only the reasonably wide windows
            wide_enough = (ga_dt >= 0.5/24.0) # [days]
            self.ga_t0 = ga_t0[wide_enough]
            self.ga_dt = ga_dt[wide_enough]

    def fill_ga_times(self, dt_generator):
        r'''Fill open GA intervals with synthetic GA observations.'''
        # all GA arrival times and lengths
        xobA_t0, xobA_dt = [], []
        for t0, dt in zip(self.ga_t0, self.ga_dt):
            # choose obs starting at t0, to fill up dt units of time
            xob_dt = []
            while sum(xob_dt) < dt:
                xob_dt.append(dt_generator())
            # adjust last obstime to fit in the box...
            xob_dt[-1] = dt - sum(xob_dt[:-1])
            # ...and if the last one is too small (time in days), forget it
            if xob_dt[-1] < 4.0/(24.0*60.0):
                xob_dt = xob_dt[:-1]
            # make arrival times to go with these dt's
            xob_t0 = [t0 + sum(xob_dt[:i]) for i in range(len(xob_dt))]
            # append this block of times to the big list
            xobA_t0.extend(xob_t0)
            xobA_dt.extend(xob_dt)
        # split the big list into one for each subsidiary instrument
        # (in an essentially arbitrary way)
        xob1_t0, xob1_dt = [], []
        xob2_t0, xob2_dt = [], []
        count = 0.0
        for t0, dt in zip(xobA_t0, xobA_dt):
            # will produce runs of length ~ 1/0.11 = 9
            count += 0.11
            if int(count) % 2 == 1:
                xob1_t0.append(t0)
                xob1_dt.append(dt)
            else:
                xob2_t0.append(t0)
                xob2_dt.append(dt)
        # save in object
        self.xob1_t0 = np.array(xob1_t0)
        self.xob1_dt = np.array(xob1_dt)
        self.xob2_t0 = np.array(xob2_t0)
        self.xob2_dt = np.array(xob2_dt)


########################################
###
###  Plot Container Class
###
########################################

class plotTimelineContainer(object):
    """A container class for methods relating to plotting DRM timelines
    """
    def __init__(self, args=None):
        self.args = args

    def style_plot(self):
        """Styles plots.
        """
        plt.rc('axes', linewidth=2)
        plt.rc('lines', linewidth=2)
        plt.rcParams['axes.linewidth'] = 2
        # plt.rc('font', weight='bold')

    def panel_legend(self, plot_attrs, sim):
        r'''Deduce a legend from what we plot.  Return the legend, and a tabularized version too.'''
        texts = []
        infos = []
        for plot_num in range(len(plot_attrs)):
            # unpack
            name, pos, width, color, _ = plot_attrs[plot_num]
            t0name, dtname = _
            # select obs within trange
            t0_all = getattr(sim, t0name)
            dt_all = getattr(sim, dtname)
            extent = np.sum(dt_all)
            proportion = 100.0 * (extent / sim.missionLife)
            texts.append('%s: %.0f days = %.1f%%' % (name.replace('-', ''), extent, proportion))
            infos.append(','.join((name.replace('-', ''), '%.0f' % extent, '%.1f' % proportion)))
        return ('\n'.join(texts), '\n'.join(infos))

    def plot_timeline_panel(self, sim, slew_debug=False, debug_out=False):
        # plot attributes for a "collection"
        plot_attributes = [
            ('Detection',     25, 8, ('blue', 'lightblue'),      ('det_t0',  'det_dt' )),
            ('Spectra',       15, 8, ('green', 'lightgreen'),    ('char_t0', 'char_dt_pad')),
            ('-Slew',         15, 2, ('lightgray', 'darkgray'),  ('slew_t0', 'slew_dt')),
            ('Other',          5, 8, ('chocolate', 'peachpuff'), ('ga_t0',   'ga_dt'  ))]

        # this makes an extra line for slews - I ended up using a different method though.
        # slew_debug_attrs = [
        #    ('',              20, 2, ('wheat', 'tan'),       ('slew_t0', 'slew_dt')),
        #    ]
        slew_debug_attrs = []
        
        if slew_debug:
            plot_attributes.extend(slew_debug_attrs)
        # one figure for whole plot
        fig = plt.figure(figsize=(18,20))
        # subplots starting at t_first days, each covering t_cover days
        t_first = 0 # [days]
        t_cover = 365 # [days]
        for n_plot in range(5):
            ax = fig.add_subplot(5, 1, n_plot+1)
            t_str = 'Year %d' % (n_plot+1)
            self.plot_timeline(ax, plot_attributes, sim, t_str,
                                   (t_first+n_plot*t_cover, t_first+(n_plot+1)*t_cover))

        plt.subplots_adjust(hspace=0.5)
        # add text
        if True:
            text, infos = self.panel_legend(plot_attributes, sim)
            plt.figtext(0.02, 0.05, text, verticalalignment='center', weight='bold', fontsize=16)
            if debug_out:
                fname = 'obs-timeline-info'
                fn = self.args.outpath % (fname, 'csv')
                with open(fn, 'w') as f:
                    f.write(infos)
        # show the plot
        plt.show(block=False)
        fname = 'obs-timeline'
        plt.savefig(self.args.outpath % (fname, 'png'))
        # plt.savefig(self.args.outpath % (fname, 'pdf'))

        plt.close()

    def plot_timeline_collection(self, sim):
        # plot attributes for a "collection"
        plot_attributes = [
            ('Coronagraph',   35, 8, ('blue', 'lightblue'),     ('det_t0',  'det_dt' )),
            ('Starshade',     25, 8, ('green', 'lightgreen'),   ('char_t0', 'char_dt')),
            ('-Slew',         25, 2, ('lightgray', 'darkgray'), ('slew_t0', 'slew_dt')),
            ('HWC Parallel',  12, 2, ('khaki'),                 ('all_t0',  'all_dt' )),
            ('HWC Dedicated', 15, 8, ('gold', 'goldenrod'),     ('xob1_t0', 'xob1_dt')),
            ('UVS Parallel',   2, 2, ('thistle'),               ('all_t0',  'all_dt' )),
            ('UVS Dedicated',  5, 8, ('purple', 'violet'),      ('xob2_t0', 'xob2_dt'))]
        # plots starting at t_first days, each covering t_cover days
        t_first = 0 # [days]
        t_cover = 100 # [days]
        for n_plot in range(19):
            # new figure for each plot
            fig = plt.figure(figsize=(20,5))
            ax = fig.add_subplot(111)
            t_str = '%02d' % n_plot
            self.plot_timeline(ax, plot_attributes, sim, 'Segment ' + t_str,
                                   (t_first+n_plot*t_cover, t_first+(n_plot+1)*t_cover))
            # show the plot
            plt.show(block=False)
            fname = 'day-in-the-life-seg-' + t_str
            plt.savefig(self.args.outpath % (fname, 'png'))
            # plt.savefig(self.args.outpath % (fname, 'pdf'))
            plt.close()

    def plot_timeline(self, ax, plot_attrs, sim, lbl, t_range):
        # plot parameters
        #   title, y-pos, y-width, color(s), (t0_attr, dt_attr)
        # put up the plots
        for plot_num in range(len(plot_attrs)):
            # unpack
            name, pos, width, color, _ = plot_attrs[plot_num]
            t0name, dtname = _
            # select obs within trange
            t0_all = getattr(sim, t0name)
            dt_all = getattr(sim, dtname)
            t_index = np.logical_and((t0_all+dt_all) >= t_range[0], t0_all <= t_range[1])
            self.plot_timeline_single(ax, pos, width, color, t0_all[t_index], dt_all[t_index])

        ## plot styling
        self.style_plot()
        # plot labeling info
        y_lbl = [attrs[0] for attrs in plot_attrs if not attrs[0].startswith('-')]
        y_pos = [attrs[1] for attrs in plot_attrs if not attrs[0].startswith('-')]
        y_del = max([attrs[2] for attrs in plot_attrs if not attrs[0].startswith('-')])/2.0 # half-height
        # plot labeling
        ax.set_xlim(t_range)
        ax.set_ylim((min(y_pos)-y_del-1, max(y_pos)+y_del+1)) # pad by 1
        ax.xaxis.set_tick_params(labelsize=16)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(y_lbl, fontsize=16)
        ax.set_xlabel('Time Since Mission Start [day]', weight='bold', fontsize=16)
        if sim.show_seed:
            title = 'Mission Timeline for %s: %s' % (sim.name, lbl)
        else:
            title = 'Mission Timeline: %s' % (lbl, )
        plt.title(title, weight='bold',fontsize=18)

    def plot_timeline_single(self, ax, center, width, color, t0, dt):
        r'''Helper function to plot one set of bars for a single instrument.'''
        # re-order as: [(t0, dt), (t0, dt), ... (t0, dt)]
        time_arg = np.transpose(np.stack((t0, dt))) 
        ax.broken_barh(time_arg, (center-width/2, width), facecolors=color)


###
### possibly-useful stuff relating to Exosims time parameterization
###

# allModes = outspec['observingModes']
# mode1 = [mode for mode in allModes if 'detectionMode' in mode.keys() or 'detection' in mode.keys()]
# assert len(mode1) >= 1, 'This needs to be enhanced'
# mode = mode1[0]
# if not 'timeMultiplier' in mode.keys():
#     mode['timeMultiplier'] = 1.

# arrival_times = [drm[i]['arrival_time'] for i in LD]
# sumOHTIME = outspec['settlingTime'] + outspec['starlightSuppressionSystems'][0]['ohTime']
# det_times = [drm[i]['det_time'].value*(mode['timeMultiplier'])+sumOHTIME for i in LD]
# det_timesROUNDED = [round(drm[i]['det_time'].value*(mode['timeMultiplier'])+sumOHTIME,1) for i in LD]
# ObsNums = [drm[i]['ObsNum'] for i in LD]
# char_times = [drm[i]['char_time'].value*(1.+outspec['charMargin'])+sumOHTIME*(drm[i]['char_time'].value > 0.) for i in LD]


def main(args):
    # load information for a generator of plausible observation times
    if args.dayinthelife:
        values, probs = load_dt_generator_info(args.obstime_csv)
        dt_generator = lambda: draw_obstime(values, probs)
    # load the simulation run
    sim = SimulationRun(args.drm, args.script)
    sim.get_obs_times()
    # find GA time windows
    sim.find_ga_times()
    # fill GA time windows with other observations
    if args.dayinthelife:
        sim.fill_ga_times(dt_generator)
    # open a plot container object
    plotter = plotTimelineContainer(args)
    # make plots
    if args.dayinthelife:
        plotter.plot_timeline_collection(sim)
    if args.synoptic:
        plotter.plot_timeline_panel(sim, slew_debug=args.slew_debug, debug_out=args.debug_out)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Plot timeline of observations reported in DRMs.", epilog='')
    parser.add_argument('script', metavar='SCRIPT', default='', help='json script file')
    parser.add_argument('drm', metavar='DRM', default='', help='drm file')
    parser.add_argument('-t', default='', type=str,
                            dest='obstime_csv', help='GA duration CSV file. Default = builtin values. Un-needed if -d given.')
    parser.add_argument('-o', default='./%s.%s', type=str, dest='outpath', help='Output file pattern.')
    #parser.add_argument('-v', default=False, action='store_true', 
    #                        dest='verbose', help='Verbosity.')
    parser.add_argument('-s', default=True, action='store_false', 
                            dest='synoptic', help='Omit synoptic plot.')
    parser.add_argument('-d', default=True, action='store_false', 
                            dest='dayinthelife', help='Omit day-in-the-life plot series.')
    parser.add_argument('--slew_debug', default=False, action='store_true', 
                            help='Show extra information to debug slews.')
    parser.add_argument('--debug_out', default=True, action='store_false', 
                            help='Suppress extra information usually written to a CSV file.')
    args = parser.parse_args()
    args.progname = os.path.basename(sys.argv[0])
    
    if args.outpath.count('%s') != 2:
        sys.stderr.write('Fatal.  Outpath (%s) must have two %%s patterns.' % args.outpath)
        sys.exit(1)

    # set umask in hopes that files will be group-writable
    os.umask(0o002)

    # ensure enclosing dir exists
    directory = os.path.dirname(args.outpath % ('dummy', 'txt'))
    tries = 0
    while not os.path.exists(directory) and tries < 1:
        try:
            os.makedirs(directory)
        except OSError:
            tries += 1
            time.sleep(random.uniform(0.1, 0.2))
            pass # hope it was just concurrent creation 

    # seed the RNG, which is used for generating GA-instrument obs times
    # we use a deterministic seed for repeatability
    np.random.seed(12345)

    # do it
    main(args)
    sys.exit(0)

