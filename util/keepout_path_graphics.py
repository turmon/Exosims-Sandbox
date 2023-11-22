#!/usr/bin/env python
#
# keepout_path_graphics: Plot targets, planets, and keepout regions, with optional DRM overlay
#
# For usage, use the -h option.  Some options may be described there but not here.
#
# Produces a time-stepped map, in sky coordinates, of stellar targets, planets, keepout regions,
# and, if a DRM is provided (see below), detections, characterizations, and slews.
# Output is to a movie (.mp4), or individual frames (.png).  Also, some cumulative keepout
# summary statistics, and plots, can optionally be generated.
#
# DRM or observing-tour display:  If pickles summarizing a DRM are supplied (--drm), then the
# observing tour is loaded and the observations are shown. It is assumed that --drm and --spc
# are used together. 
#
# Typical usage:
#   keepout_path_graphics.py -s 0 -l 0.2 -d 0.5 -m $HOME/keepout.mp4 ./sampleScript_coron.json 
# where:
#  -s is the start-time of the movie
#        if zero, the start-time from the script is used.
#        if a float in [0,1], the movie is started that proportion of the way through the mission.
#        if a number <= 10000 is given, this is the offset from mission start in days
#        if a number > 10000 is given, this is the MJD start time.
#  -l is the length in years
#        if zero is given (-l 0), the duration from the script is used.
#  -d is the delta-t between frames in days
#  -m is the name of the movie
# optionally:
#  -e     -- use equatorial coordinates as opposed to ra/dec, HIGHLY recommended
#  -f DIR -- the directory name for .png frame-by-frame output
#  -c DIR -- the directory name for cumulative keepout output
#  --drm FILE -- the file containing the DRM as a pickle
#  --spc FILE -- the file containing the SPC as a pickle
#
# experts only:
#  -x SCRIPT -- the given SCRIPT filename is loaded on top of the argument script
#               if the given SCRIPT name begins with !, it is treated as a
#               json literal rather than a filename


# history:
#  Christian Delacroix, Oct 2016
#  Michael Turmon, Oct 2016
#  Christian Delacroix, Mar 2017
#  Christian Delacroix, May 2017 (added occulter koAngleMax)
#  turmon may 2017: plot cumulative results
#  turmon aug 2017: add DRM overlay
#  turmon nov 2018: add occulter keepout to existing detector keepout
#  turmon mar 2020: python3, astropy4
#  turmon dec 2020: updated to handle starshade-only cases better
#  turmon 2023: simplify, fix multiple issues with newer DRM formats
#
# Note: The "cumulative observability" map-format plots are fragile, because
# the "koMap" returned by Obs.keepout() changed after this code was
# written.  It formerly contained 1 map (Stars X Times), but now contains
# NS such maps, where NS is the number of StarlightSuppressionSystems.  
# Sometimes NS=1 (coronagraph-only, or starshade-only) and sometimes NS=2
# (coronagraph + starshade).
# The code here isn't as smart as it could be about this.  It assumes 
# if there is a starshade, then koMap[-1,:,:] has its keepout.  Our starshades 
# happen to be the last observing mode, so this works, but it's fragile.
# I think this would be fixed by looking at observingModes to find the right index.
#
# TODO: This was the first keepout-extracting code we wrote.  It does not separate
# compilation-of-information from generation-of-graphics.  So it has a main loop 
# that is lengthy and unfactored (make_graphics(), about 500 lines)


import sys
import glob
import argparse
import os
import os.path
import json
import six.moves.cPickle as pickle
import random
import time
from collections import Counter
import EXOSIMS
import EXOSIMS.MissionSim
import numpy as np
import astropy
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord
import matplotlib; matplotlib.use('Agg') # not interactive: don't use X backend
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as patches
import matplotlib.lines as lines
from matplotlib import collections as mplc
from matplotlib import rcParams

# For tqdm.__version__ >= 4.66, you can:
#   set environment variable TQDM_DISABLE=1
#   os.environ["TQDM_DISABLE"] = "1"
# 
# To allow older versions, and to keep this self-contained,
# I'm patching tqdm (used in Obs.keepout) like this:
from tqdm import tqdm
from functools import partialmethod
tqdm.__init__ = partialmethod(tqdm.__init__, disable=True)


# keeps xlabel from being chopped off
rcParams.update({'figure.autolayout': True})

# graphics styles
axlabel_style = dict(fontsize=10, fontweight='bold')
title_style = dict(fontsize=12, fontweight='bold')
tick_style = dict(labelsize=10)

# unpickling python2/numpy pickles within python3 requires this
PICKLE_ARGS = {} if sys.version_info.major < 3 else {'encoding': 'latin1'}

############################################################
#
# Utility Functions
#
############################################################

def ensure_dir(directory):
    r'''Ensure enclosing dir exists, especially in multiprocessing context.'''
    # the loop guards against another process creating the dir after we checked,
    # but before the os.makedirs() call.
    tries, done = 0, False
    while not os.path.exists(directory) and not done:
        try:
            os.makedirs(directory, 0o775)
            done = True # ensure we exit the while
        except OSError:
            # hope it was just concurrent creation
            tries += 1
            time.sleep(random.uniform(0.1, 0.2))
            if tries == 2:
                raise # give up

def get_koangles(OS):
    r'''Get koangles array from information in the OpticalSystem.
    The (un-pythonic) code here is copied from Prototypes/SurveySimulation so that 
    we generate exactly the koangles as expected by the caching mechanism and
    the Observatory.keepout() method.'''
    # choose observing modes selected for detection (default marked with a flag)
    allModes = OS.observingModes

    nSystems  = len(allModes)
    systNames = np.unique([allModes[x]['syst']['name'] for x in np.arange(nSystems)]).tolist()
    koStr     = ["koAngles_Sun", "koAngles_Moon", "koAngles_Earth", "koAngles_Small"]
    koangles  = np.zeros([len(systNames),4,2])
    tmpNames  = list(systNames)
    cnt = 0
    for x in np.arange(nSystems):
        name = allModes[x]['syst']['name']
        if name in tmpNames:
            koangles[cnt] = np.asarray([allModes[x]['syst'][k] for k in koStr])
            cnt += 1
            tmpNames.remove(name)
    return koangles

def extract_time(tm):
    r'''Unbox a time, if it is boxed, for compatibility with multiple DRMs.'''
    try:
        return tm.value
    except AttributeError:
        return tm

def get_char_status(obs):
    r'''Utility function, gets char status from a drm observation.
    '''
    if 'char_status' in obs:
        # starshade char - has a unitary char_status key
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
    
def strip_units(x):
    r'''Strip astropy units from x.'''
    # TODO: allow coercing units to a supplied value
    if hasattr(x, 'value'):
        return x.value
    else:
        return x

def get_slew_begin(obs):
    r'''Return the time of the beginning of the slew, if any, for obs'''
    slew_time = strip_units(obs.get('slew_time', 0.0))
    arrival_time = strip_units(obs['arrival_time'])
    return arrival_time - slew_time


class GraphicsStyle(object):
    r'''Singleton class that holds graphics styles.'''
    # mapping of conditions to colors
    # formerly: char_ok   = (0.0, 0.5, 1.0) # sky blue
    char_HZ   = (0.0, 0.7, 0.0) # green
    char_ok   = (0.5, 0.2, 1.0) # purple
    char_fail = (0.8, 0.1, 0.0) # brick-red
    det_HZ    = (0.0, 0.7, 0.0) # green
    det_ok    = (0.5, 0.2, 1.0) # purple
    det_fail  = (0.8, 0.1, 0.0) # brick-red
    # line/shape styles
    char_mec  = 'black' # characterization, marker edge color
    det_mec = None # detection, marker edge color
    char_lw = 0.5 # characterization, linewidth
    det_lw = 0.0 # detection, linewidth

    def __init__(self): pass

    def count2pt(self, count):
        r'''Display "points" (i.e., circle radius in point units) per detection.'''
        return np.sqrt(0.3 + 1.7 * count)

    def count2pt2(self, count):
        r'''Display "points" (i.e., markersize point^2 units) per detection.
        The relationship between count2pt vs. count2pt2 is empirical.'''
        return 1.7 + 1.0 * count


# Container class for loading canned observational data from the DRM
class ObserveInfo(object):
    r"""Observations made by a mission, in order, as loaded from an external DRM pickle."""

    def load_from_spc_file(self, loc, spc):
        r'''Load the DRM, and star/planet info from a "spc" file given as an argument.
        This drm + spc file transfer is compatible the Exosims ipyparallel output.'''
        # these will fail noisily if there is no file present
        print('Loading DRM from', loc)
        self.drm = pickle.load(open(loc, 'rb'), **PICKLE_ARGS)
        # spc file contains a dict with many fields - save them all
        self.spc = pickle.load(open(spc, 'rb'), **PICKLE_ARGS)
        # and then extract the ones we need
        self.snm = self.spc['Name']
        #self.sco = self.spc['coords'] # we have the script, so coords are not needed
        #self.spc_coord = self.sco.heliocentrictrueecliptic
        # planet Equivalent Insolation Distance (EID), for HZ test
        self.EID_planet = (self.spc['a'] / np.sqrt(self.spc['L'][self.spc['plan2star']])).value
        # planet radius, in EarthRad
        self.Rp_planet = self.spc['Rp'].value

    def __init__(self, loc, spc):

        # load DRM and Star-Planet info
        if spc:
            self.load_from_spc_file(loc, spc)
        # number of observations
        N = len(self.drm)
        # stars we visited, in DRM order
        self.drm_s_ind = [r['star_ind'] for r in self.drm]
        # visit times
        self.times = np.array([extract_time(r['arrival_time']) for r in self.drm])
        # summaries
        # Note, there is room for improvement with recording successful observations,
        # and HZ observations.  Specifically, we could record number of dets/chars/HZ,
        # or we could record only unique dets/chars/HZ.  It really depends on what we
        # want to plot.
        self.was_det  = np.zeros((N,), dtype=bool) # detection or char?
        self.success  = np.zeros((N,), dtype=bool) # was it successful for any planet?
        self.hab_zone = np.zeros((N,), dtype=bool) # were there successful HZ observations?
        for i, obs in enumerate(self.drm):
            # deal with various key possibilities
            if 'det_status0' in obs:
                Dkey = 'det_status0'
            elif 'det_status1' in obs:
                Dkey = 'det_status1'
            else:
                Dkey = 'det_status'
            # was the observation a detection?
            was_det = (Dkey in obs)
            self.was_det[i] = was_det
            # query detections/characterizations for status in the same way
            if was_det:
                obs_status = obs[Dkey]
            else:
                obs_status = get_char_status(obs)
            # 03-Nov-2021, turmon // Nov-2023, turmon, believe the below is obsolete
            # inserted due to some DRM's containing records that are neither det nor char
            ##    print('Masked error: set DRM[%d][%s] = 0 because key is absent (star_ind = %d)' %
            ##              (i, st_key, obs['star_ind']))
            # was the observation successful?
            self.success[i] = np.any(np.array(obs_status) > 0)
            # was any successful observation done on a HZ planet?
            # 1/ planet-indexes of successful observations, if any
            success_plan_inds = np.array(obs['plan_inds'],dtype=int)[np.where(obs_status > 0)[0]]
            # 2/ indicator for each successful planet, Is the planet in the HZ?
            #    turmon 11/2021: changed the upper boundary of L_planet from 1.55 to 1.67
            hab_zone = np.logical_and(self.EID_planet[success_plan_inds] >= 0.40,
                                      self.EID_planet[success_plan_inds] <= 1.67)
            # 3/ was any successful planet observation in the HZ?
            #    Note: perhaps we actually want a count?
            self.hab_zone[i] = np.any(hab_zone)

    def tour_summary(self, args):
        '''Compute and optionally dump tour summary statistics.'''
        Nstar = len(self.snm)
        tour = [obs for i, obs in enumerate(self.drm) if not self.was_det[i]]
        # list of all stars visited
        self.visited = np.array([d['star_ind'] for d in tour], dtype=int)
        # total cumulative number-of-visits, by star
        self.vcount = np.bincount(self.visited, minlength=Nstar)
        # slew adjacency matrix
        self.visit2 = np.zeros((Nstar, Nstar), dtype=int)
        # slew transition list (on what chars did we go from d1->d2)
        visit_when = np.empty((Nstar, Nstar), dtype=object)
        for d1 in range(Nstar):
            for d2 in range(Nstar):
                visit_when[d1,d2] = []
        # if only one thing in tour, there was no slew
        if len(tour) > 1:
            for when, d1, d2 in zip(list(range(len(tour)-1)), tour[0:-1], tour[1:]):
                self.visit2[d1['star_ind'], d2['star_ind']] += 1
                visit_when[ d1['star_ind'], d2['star_ind']].append(when)
    
        # dump some info - should use a csv-writer
        if args.out_cume:
            # visit counts
            fn = os.path.join(args.out_cume, 'path-visits.csv')
            with open(fn, 'w') as f:
                arr = ['name', 'visit_mean', 'lon', 'lat', 'dist']
                f.write(','.join(arr) + '\n')
                for i, ct in enumerate(self.vcount):
                    arr = [
                        '%s' % self.snm[i],
                        '%.1f' % ct,
                        '%.2f' % self.xpos[i],  # [deg]
                        '%.2f' % self.ypos[i],  # [deg]
                        '%.2f' % self.distance[i], # [pc]
                        ]
                    f.write(','.join(arr) + '\n')
                print("Visits written to `%s'" % fn)
            # slew counts
            fn = os.path.join(args.out_cume, 'path-slews.csv')
            with open(fn, 'w') as f:
                pos1, pos2 = np.where(self.visit2 > 0)
                arr = ['source', 'dest', 'slews', 'label']
                f.write(','.join(arr) + '\n')
                for p1, p2 in zip(pos1, pos2):
                    # e.g., "Slew Numbers: 17; 39"
                    label_string = ('"Slew Number%s: %s"' % (
                        '' if len(visit_when[p1,p2]) == 1 else "s",
                        '; '.join(str(when+1) for when in visit_when[p1,p2])))
                    arr = [
                        '%d' % p1,
                        '%d' % p2,
                        '%d' % self.visit2[p1,p2],
                        label_string
                        ]
                    f.write(','.join(arr) + '\n')
                print("Slews written to `%s'" % fn)

    def summary(self):
        s = f'Loaded DRM:\n  {len(self.drm)} observations\n  {len(set(self.drm_s_ind))} stars observed\n  {len(self.snm)} stars total'
        return s

    def lookup(self, day):
        r"""Return index of observation closest to day."""
        return np.argmin(np.abs(self.times - day))

    def obs_history(self, day):
        r"""Return observation history up to the given day: counts of detections and charaterizations.
        Each observation fits into one of the cases below.
        Coded as below, where each is a Counter indexed by star number. 
          dets0 = observation with no detections.
          dets1 = observation with >=1 detection.
          chars0 = characterization attempt, failed.
          chars1 = characterization attempt, success.
        Could enhance to add another counter, for habitable-zone detections vs. non-HZ detections.
        """
        index = self.lookup(day)
        # TODO: prefer a single dict of Counter(), instead of multiple Counter's
        detsHZ  = Counter(); dets1  = Counter(); dets0  = Counter()
        charsHZ = Counter(); chars1 = Counter(); chars0 = Counter()
        for i in range(index):
            # star index for this observation
            star_ind = self.drm[i]['star_ind']
            # FIXME: Some DRMs have chars and dets in the same observation
            if self.was_det[i]:
                # detection
                # FIXME: HZ detections will wipe out non-HZ detections
                # could fix this with a count of HZ and non-HZ detections
                if self.hab_zone[i]:
                    detsHZ[star_ind] += 1
                elif self.success[i]:
                    dets1[star_ind] += 1
                else:
                    dets0[star_ind] += 1
            else:
                # characterization
                if self.hab_zone[i]:
                    charsHZ[star_ind] += 1
                elif self.success[i]:
                    chars1[star_ind] += 1
                else:
                    chars0[star_ind] += 1
        return detsHZ, dets1, dets0, charsHZ, chars1, chars0

    def obs_summary(self, day, normal, z0):
        r"""Return a summary of observation history up to given day.

        If normal is False, it is the last day, and slightly different summary rules hold.
        The summary is a list of tuples, intended to be: star_number, color, size, extra.
        In this case, `extra' is a dictionary of point properties, which will be given 
        to the matplotlib circle-plotter.
        To control plot order, we need to have the base z-order (z0) of the returned
        summaries.  When we need objects to come out on top, we use z0+(dz), as follows:
          dets: z0; dets-on-top: z0+1; chars: z0+2; chars-on-top: z0:3
        Because of closely-packed binaries, it is important to put chars on top.
        """
        # observations to emphasize:
        #  set of star indexes in the present slew
        stars_important = set(self.stars_on_slew(day))
        # obtain graphics styles
        GS = GraphicsStyle()

        # characterization marker style - always on top of dets
        c_style = dict(ec=GS.char_mec, lw=GS.char_lw, zorder=z0+2)
        # make rv, the list of observations-to-date
        rv = []
        detsHZ, dets1, dets0, charsHZ, chars1, chars0 = self.obs_history(day)
        for star in (detsHZ + dets1 + dets0 + charsHZ + chars1 + chars0).keys():
            # if not "important", detected stars will be dimmer
            # if last frame (not normal), no detected stars will be dimmed
            # in either case, set no edge: the edge overlaps the body of the circle,
            # and then, when made transparent, donut shapes result.
            if (star in stars_important) or (not normal):
                d_style = dict(ec=GS.det_mec, lw=GS.det_lw, zorder=z0)
            else:
                d_style = dict(ec=GS.det_mec, lw=GS.det_lw, zorder=z0, alpha=0.3)
            # order of *if* statements encodes display preference order - show only one glyph per star
            if chars1[star] > 0 and charsHZ[star] > 0:
                # HZ and regular: circle within a circle. Omit fails.
                # count2pt(1+...) gives extra room for the bounding circle
                rv.append( (star, GS.char_HZ, GS.count2pt(1+chars1[star]+charsHZ[star]), c_style) )
                rv.append( (star, GS.char_ok, GS.count2pt(chars1[star]), dict(c_style, zorder=z0+3)) )
            elif charsHZ[star] > 0:
                rv.append( (star, GS.char_HZ, GS.count2pt(charsHZ[star]), c_style) )
            elif chars1[star] > 0:
                rv.append( (star, GS.char_ok, GS.count2pt(chars1[star]),  c_style) )
            elif chars0[star] > 0:
                rv.append( (star, GS.char_fail, GS.count2pt(chars0[star]), c_style) )
            elif dets1[star] > 0 and detsHZ[star] > 0:
                # HZ and regular: circle within a circle. Omit fails.
                rv.append( (star, GS.det_HZ, GS.count2pt(1+dets1[star]+detsHZ[star]), d_style) )
                rv.append( (star, GS.det_ok, GS.count2pt(dets1[star]), dict(d_style, zorder=z0+1)) )
            elif detsHZ[star] > 0:
                rv.append( (star, GS.det_HZ, GS.count2pt(detsHZ[star]), d_style) )
            elif dets1[star] > 0:
                rv.append( (star, GS.det_ok, GS.count2pt(dets1[star]), d_style) )
            elif dets0[star] > 0:
                # all detections failed
                rv.append( (star, GS.det_fail, GS.count2pt(dets0[star]), d_style) )
            else:
                assert False, 'Condition should be unreachable'
        return rv

    def char_bracket(self, day):
        r"""Return pair of characterization observations bracketing a given day."""
        center = self.lookup(day)
        hi = np.min(np.where(self.was_det[center:] == False)[0])
        # char_status is not at top-level of all the DRM's
        # lo = max([i for i in range(center) if 'char_status' in self.drm[i]])
        lo = max([i for i in range(center) if not self.was_det[i]])
        if not hi or not lo:
            return None, None
        else:
            return lo, hi
            
    def char_before(self, day):
        r'''Return index of the previous characterization BEFORE day'''
        center = self.lookup(day)
        try:
            endpoint = np.where(self.was_det[0:center] == False)[0][-1]
        except IndexError:
            endpoint = 0
        return endpoint

    def char_after(self, day):
        r'''Return index of the next characterization AFTER day'''
        center = self.lookup(day)
        try:
            endpoint = center + np.where(self.was_det[center:] == False)[0][0]
        except IndexError:
            endpoint = len(self.was_det) - 1
        return endpoint

    def char_to_date(self, day):
        r"""Return all characterization observations up to a given day."""
        # catch the next characterization AFTER day
        endpoint = self.char_after(day)
        # return DRM entries from 0:endpoint-1 that were NOT detections
        return np.take(self.drm, np.where(self.was_det[:endpoint+1] == False)[0])
        
    def stars_on_slew(self, day):
        r"""Return list of all stars in the characterization-slew-interval containing day."""
        i0 = self.char_before(day)
        i1 = self.char_after(day)
        return self.drm_s_ind[i0:i1+1]


def lonlat2xyz(lon, lat):
  x = np.cos(np.radians(lat)) * np.cos(np.radians(lon))
  y = np.cos(np.radians(lat)) * np.sin(np.radians(lon))
  z = np.sin(np.radians(lat))
  return (x,y,z)

def xyz2lonlat(x, y, z):
  r = np.sqrt(x*x + y*y + z*z)
  lon = np.mod(np.degrees(np.arctan2(y, x)), 360)
  lat = np.degrees(np.arcsin(z / r))
  return lon, lat

def gcp(lon1, lat1, lon2, lat2):
    r"""GCP = great circle points from (lon1,lat1) -> (lon2,lat2).

    Point-count is auto-determined for smoothness without waste."""
    Npt = 200
    p1 = np.array(lonlat2xyz(lon1, lat1))
    p2 = np.array(lonlat2xyz(lon2, lat2))
    delta = p2 - p1
    # length of the line: between 0 and 2, generally < 1
    r = np.linalg.norm(delta)
    d = np.linspace(0, 1.0, int(Npt*r + 10))
    # sample points along the line connecting p1 and p2
    pts = p1 + np.outer(d, delta)
    # convert these points to lon,lat
    pts_ll = [xyz2lonlat(p[0], p[1], p[2]) for p in pts]
    return pts_ll

def wrap_paths(p1, p2):
    r'''Construct a path from p1 to p2, allowing for wrap-around of longitude.

    Note that latitude does not wrap, because +90 is distinct from -90.'''
    # destructure the pairs
    x1, y1 = p1
    x2, y2 = p2

    # Construct a great-circle path from p1 to p2: a series of segments
    #   for reference, the one-segment path from p1 -> p2 is:
    #     lines = [[(x1, y1), (x2, y2)]]
    pts = gcp(x1, y1, x2, y2)
    lines = []
    # allow for longitude wrapping
    for i in range(len(pts)-1):
        p1 = pts[i]
        p2 = pts[i+1]
        # bisect line segments that wrap around
        if abs(p1[0] - p2[0]) > 90:
            # a midway latitude, not distance-weighted: ok for short segment
            midlat = (p1[1] + p2[1]) * 0.5
            if p1[0] > p2[0]:
                # 360 -> 0 wrap
                lines.append([(p1[0], p1[1] ), (360.0, midlat)])
                lines.append([(  0.0, midlat), (p2[0], p2[1] )])
            else:
                # 0 -> 360 wrap
                lines.append([(p1[0],  p1[1]), (  0.0, midlat)])
                lines.append([(360.0, midlat), (p2[0], p2[1] )])
        else:
            lines.append([(p1[0], p1[1]), (p2[0], p2[1])])
    return lines

def get_arrow(lines, length):
    r""""Pick one of the line segments within lines, and make an arrow point along it."""
    n0 = len(lines) - 1
    arrow_0 = np.array(lines[n0][0])
    arrow_1 = np.array(lines[n0][1])
    delta = arrow_1 - arrow_0
    # make it unit norm - safeguard in case line is degenerate
    delta = delta / np.maximum(1e-6, np.linalg.norm(delta))
    return arrow_0, length * delta

def wrap_path(x, y):
    r'''Given a path as lists of x's and y's, separate it into contiguous path segments.

    The input paths are lists of (lon,lat), and if the jump between adjacent values is
    too large, this indicates a jump such as a wrap-around from 359.9 -> 0.1 in longitude.
    Think of a circle in (lon,lat) near the right edge of the plot.
    We separate the path at these points.  The return value looks like:
       lines = [[(x1, y1), (x2, y2), ...], [(xP, yP), (xP+1, yP+1), ...] , ... ]
    We also return a boolean, "closed", indicating if the resulting line segments
    each enclose a closed boundary.'''
    lines = []
    line1 = []
    for i in range(len(x)):
        line1.append((x[i], y[i]))
        # index of next point, wrapping around
        i_next = i+1 if i+1 < len(x) else 0
        # start a new line segment if the next point will cause a big jump
        # note, we can jump in longitude, but not in latitude
        if abs(x[i] - x[i_next]) > 90:
            # finish the line appropriately depending on how the wrapping goes
            midpoint = (y[i] + y[i_next]) * 0.5
            if x[i] > 340 and x[i_next] < 20:
                # wrapping upward from 360 -> 0: cut line at 360, start at 0
                line1.append((360.0, midpoint))
                lines.append(line1)
                line1 = [(0.0, midpoint)]
            elif x[i] < 20 and x[i_next] > 340:
                # wrapping down from 0 to 360: cut at 0, start at 360
                line1.append((0.0, midpoint))
                lines.append(line1)
                line1 = [(360.0, midpoint)]
            else:
                # harder cases: near poles, longitude can jump but not wrap,
                # so in this case make no assumption: just delete the jump
                lines.append(line1)
                line1 = []
    # leftovers
    if len(line1) > 0:
        lines.append(line1)
    # perform some ad hoc shape joins
    if len(lines) == 3:
        # this special case joins circles that have been split
        # beginning <jump> middle <jump> end
        #   ==>
        # end + beginning <jump> middle
        lines = [lines[2] + lines[0], lines[1]]
        # typically this is a circle split across lon=360, so it is two closed shapes now
        closed = True
    elif len(lines) == 2:
        lines = [lines[1] + lines[0]]
        # typically this is a circle enclosing the pole, and it is not one closed shape
        # (note, it would be possible to close this by including the polar strip)
        closed = False
    else:
        # no jumps => boundary is closed
        closed = True
    return lines, closed

def spherical_cap(lon1, lat1, theta):
    r'''Return lon-lat point-list of a spherical cap of theta degrees around (lon,lat).

    Envision the (lon,lat) point with a line connecting it to the origin of a unit sphere,
    and then a plane orthogonal to that line.  To generate the cap, we construct a basis 
    spanning the plane, and then draw a circle of radius tan(theta) in that plane.  Then
    we extract the lon/lat values of that circle, and turn them into a list of connected
    segments.  This is needed because the lon/lat circle can overflow the edge of the
    lon/lat plot, and appear as two half-circles.'''
    p_xyz = np.array(lonlat2xyz(lon1, lat1)).reshape((3,1))
    # orthogonal rotation matrix with p as one basis
    # the first column of R is p_xyz, and the other two (the ones we want) span the
    # orthongonal complement of p_xyz.  We use the latter two as a basis for the circle
    R = np.linalg.svd(p_xyz)[0]
    # radius of circle
    delta = np.tan(np.radians(theta))
    # sample [0,2pi], larger circles => more samples, very high-latitude => even more 
    twopi = np.linspace(0.0, 2*np.pi, int(theta*5+12) + (200 if theta > 50 else 0))
    # circle in the y-z plane of radius delta, with x == 0
    circle = np.stack((0*twopi, delta*np.cos(twopi), delta*np.sin(twopi)))
    # q = p + R * circle
    # rotate the circle to live in the plane perpendicular to p_xyz, and
    # translate the rotated circle to be centered about p_xyz
    q_xyz = p_xyz + np.dot(R, circle)
    # find the lon/lat corresponding to q
    lo, la = xyz2lonlat(q_xyz[0,:], q_xyz[1,:], q_xyz[2,:])
    # cut the circle if, e.g., lon jumps from 359.9 to 0.1,
    # making the circle into a list of line segments
    lola, closed = wrap_path(lo, la)
    return lola, closed


def run_max(x):
    r"""Return the maximal run of True along each row of a 2d matrix x.

    It is assumed that all rows start and end False, with one or more
    transitions to True along the way.  You can False-pad x on left
    and right to ensure this."""
    # 0->1 and 1->0 transitions, and their row,col locations in x
    row,up   = np.where(np.diff(x * 1) > 0)
    row,down = np.where(np.diff(x * 1) < 0)
    # run-lengths for each run (0, 1, or >1 runs per row of x)
    run_lengths = down - up
    # tabulate maximal run per row of x
    x_runs = np.zeros(x.shape[0])
    for r in np.unique(row):
        x_runs[r] = np.max(run_lengths[row == r])
    return x_runs
        
def circle_marker(**given_props):
    r'''Return a Line2D with circular markers for an axis legend, given_props over-rides.
    Some key properties to pass in: markeredgecolor (mec), markeredgewidth (mew), 
    markerfacecolor (mfc), markersize.'''
    # x, y are not used when in the legend
    # lw=0 => turn off lines between markers
    # markers set to circles
    # mew=0 => turn off marker edges by default
    props = dict(xdata=[0], ydata=[0], lw=0, marker='o', markeredgewidth=0)
    props.update(given_props)
    return plt.Line2D(**props)



############################################################
#
# Main routine
#
############################################################


def make_graphics(args, xspecs):

    # instantiate sim object, allowing for extra specs
    #   we use the contents of "sim" to control movie-making

    # Formerly (7/2023), th was just:
    #sim = EXOSIMS.MissionSim.MissionSim(args.script, **xspecs)
    # Now: simplify things to eliminate the dependence on
    # the JPL run_sim in Local/. So, load specs file and simplify
    # its SurveyEnsemble before object instantiation.
    try:
        with open(args.script, 'rb') as g:
            specs = json.load(g)
    except:
        sys.stderr.write('Failed to open script file "%s"\n' % args.script)
        raise
    # replace SurveyEnsemble with the Prototype
    # this mostly eliminates the need for Local/ to be import'able
    specs['modules']['SurveyEnsemble'] = " "
    # FIXME: for now (8/2023), this line removes the xspecs, which was aimed ONLY 
    # to chillax the JPL run_sim.
    # Since then, we have also removed the caller's (util/drm-to-movie.sh)
    # line giving these xspecs:
    #   -x '!{"ensemble_mode":"init-only"}'
    specs_pool = {**specs, **xspecs}

    sim = EXOSIMS.MissionSim.MissionSim(**specs_pool)
    
    ########################################
    #
    # set up parameters based on input arguments
    #
    ########################################
    # TODO: do this somewhere else

    # movie output: False, or a path (file)
    out_movie = args.out_movie
    # frame output: False, or a path (*directory*, not file)
    out_frames = args.out_frames
    # cumulative output: False, or a path (*directory*, not file)
    out_cume = args.out_cume
    # position output: False, or a path
    out_position = args.out_position
    # time-between-frames
    dt = args.delta_t
    # start time
    #   TODO: starting the movie in mid-DRM needs more to work correctly
    if args.t0 > 10000:
        # absolute MJD
        startTime = args.t0
    elif args.t0 > 0.0 and args.t0 <= 1.0:
        # fractional offset, round down to nearest day
        startTime = np.fix(sim.TimeKeeping.missionStart.value + args.t0 * (365.25 * sim.TimeKeeping.missionLife.value))
    else:
        # offset in days
        startTime = sim.TimeKeeping.missionStart.value + args.t0
    # duration
    if args.years > 0:
        nyears = args.years
    else:
        nyears = sim.TimeKeeping.missionLife.value # [years]
    # overall mission start time (vs. movie start time)
    startTimeMission = sim.TimeKeeping.missionStart.value
        
    # times are in MJD
    endTime = min(startTime + nyears*365.0, sim.TimeKeeping.missionFinishAbs.value)
    if endTime == startTime:
        # if nyears == 0, let there be one frame anyway
        currentTimes = Time(np.array([startTime]), format='mjd', scale='tai')
    else:
        currentTimes = Time(np.arange(startTime, endTime, dt), format='mjd', scale='tai')
    
    # plot movie in (equatorial) ra/dec coordinates vs. (ecliptic) lat/lon coordinates
    ra_dec_coord = args.ra_dec_coord
    if ra_dec_coord:
        # attribute name for SkyCoord's
        xattr_name, yattr_name = 'ra', 'dec'
        # text labels
        movie_xlabel, movie_ylabel = 'RA', 'DEC'
    else:
        # attribute name for SkyCoord's
        xattr_name, yattr_name = 'lon', 'lat'
        # text labels
        movie_xlabel, movie_ylabel = 'Longitude', 'Latitude'
    
    ########################################
    #
    # Filenames
    #
    ########################################
    
    # Graphics -- need to associate the movie with the figure now
    plt.close('all')
    fig = plt.gcf()
    # make aspect ratio close to 2 wide x 1 tall for a lat/lon plot
    fig.set_size_inches(6.4,3.5)
    
    if out_movie:
        # movie hooked to "fig"
        #plt.rcParams['animation.ffmpeg_path'] = '/usr/local/bin/ffmpeg'
        #FFMpegWriter = animation.writers['ffmpeg']
        # turmon 12/2021: this is fragile: path to ffmpeg on the sec383 systems
        plt.rcParams['animation.ffmpeg_path'] = '/usr/local/anaconda3/envs/python39/bin/ffmpeg'
        FFMpegWriter = animation.FFMpegWriter
        movie = FFMpegWriter(fps=15, bitrate=2500)
        movie.setup(fig, out_movie, 200) # last arg is dpi
    if out_frames:
        ensure_dir(out_frames)
    if out_cume:
        ensure_dir(out_cume)
    if out_position:
        ensure_dir(out_position)
        try:
            fp_position = open(out_position + '/position.csv', 'w')
        except IOError:
            print("Error: could not write position-file within '%s'" % out_position)
            raise
    else:
        fp_position = None
    

    ########################################
    #
    # Set up time-by-time loop
    # 
    ########################################
    
    # DRM ("Observe Info"), which can be None if it was not specified
    OI = args.OI
    # get objects
    OS = sim.OpticalSystem
    Obs = sim.Observatory
    TL = sim.TargetList
    sInds = np.arange(TL.nStars)
    evergood_coro  = np.zeros(TL.nStars) # coronagraph only
    evergood_shade = np.zeros(TL.nStars) # with starshade
    # koangles is Sx4x2: #systems x {sun, earth, moon, small bodies} x {min,max}
    # this value is used by the Obs.keepout() method
    koangles = get_koangles(OS)

    # set star coordinates from TL, rather than from the SPC file
    # FIXME: sometimes the OI has coords, and sometimes the TL has coords.
    if OI and not hasattr(OI, 'coords'):
        if False:
            # ra/dec
            # -- this is an old code path, not used any more
            OI.xpos = TL.coords.ra.value
            OI.ypos = TL.coords.dec.value
        else:
            # heliocentric lon/lat
            OI.xpos = TL.coords.barycentrictrueecliptic.lon.value
            OI.ypos = TL.coords.barycentrictrueecliptic.lat.value
            OI.distance = TL.coords.barycentrictrueecliptic.distance.value

    # Provide diagnostic summary on what modes we found
    #  coro: the first observingMode that is a detection mode
    detMode  = [mode for mode in OS.observingModes if mode['detectionMode'] == True][0]
    print('Coronagraph keepout from: %s + %s' % (detMode['systName'], detMode['instName']))
    #  shade: the first observingMode that has an occulter, or None
    try:
        shadeMode = [mode for mode in OS.observingModes if mode['syst']['occulter'] == True][0]
        print('Occulter keepout from: %s + %s' % (shadeMode['systName'], shadeMode['instName']))
    except IndexError:
        shadeMode = None
        print('No occulter observingMode, so no occulter keepout will be generated.')

    # for progress indication
    system_t0 = time.time()
    system_dt = 10.0 # seconds
    
    ########################################
    #
    # Loop over times
    # 
    ########################################
    
    # array of all kogood vectors: re-initialize within loop when we know the size
    kogoods_coro, kogoods_shade = None, None

    #if shadeMode:
    #    import pdb; pdb.set_trace()

    print('Beginning simulation propagation through time...')
    Ntime = len(currentTimes)
    for i in range(Ntime):
        # print progress if needed
        if time.time() > system_t0 + system_dt:
            sys.stdout.write('\n[%3d/%3d]' % (i+1, Ntime))
            system_t0 = time.time()
        else:
            sys.stdout.write('.')
            sys.stdout.flush()

        # final frame is different in some ways
        final_frame = (i == Ntime-1)
        
        # Build evergood array, detections
        # note: kogood is s x Ntarg x Ntime, where s = number of systems; Ntime = 1 here
        # note: koangleArr is used below, but its arrangement has changed
        currentTime = currentTimes[i]
        kogood, r_body, r_targ, _, koangleArr = Obs.keepout(TL, sInds, currentTime, koangles, returnExtra=True)
        koangleArr = koangleArr.to('deg').value
        evergood_coro += kogood[0,:,0]
        # NB: separate call not needed in new EXOSIMS
        # FIXME: this is awkward for starshade-only cases (one system, used for both det + char)
        if shadeMode:
            # -1 picks out starshade system
            evergood_shade += kogood[-1,:,0]

        # maintain kogoods
        if kogoods_coro is None:
            # initialize to false, and pad by one at time = 0 and time = end
            # these are (#targets) x (Ntime+2)
            kogoods_coro  = np.zeros((kogood.shape[1], Ntime+2), dtype=bool)
            kogoods_shade = np.zeros((kogood.shape[1], Ntime+2), dtype=bool)
        kogoods_coro[ :,i+1] = kogood[0,:,0]
        if shadeMode:
            # -1 picks out starshade system
            kogoods_shade[:,i+1] = kogood[-1,:,0]

        # strings for later titles
        time_iso = currentTime.copy()
        time_iso.format = 'iso'
        time_iso.out_subfmt = 'date_hm' # suppress sec.ms in string output below
        currentTime.out_subfmt = 'float' # reset for numeric output now
        time_mjd = '%.1f' % currentTime.mjd
        time_day = currentTime.mjd - startTimeMission
        
        # Calculate desired (ra,dec or lon,lat) coordinates of visible targets, kept-out targets, and bright bodies
        # recall that r_targ (the star positions) and r_body (the solar-system
        # body positions) are both in HE (heliocentric equatorial) coordinates.
        if astropy.version.major < 3:
            skyc_kws = dict(representation='cartesian') # this is deprecated in astropy 3+
        else:
            skyc_kws = dict(representation_type='cartesian')
        targ = SkyCoord(r_targ[0,:,0], r_targ[0,:,1], r_targ[0,:,2], **skyc_kws)
        body = SkyCoord(r_body[:,0,0], r_body[:,0,1], r_body[:,0,2], **skyc_kws)

        if ra_dec_coord:
            # SkyCoord frame to extract ra/dec
            targ = targ.heliocentrictrueecliptic.icrs
            body = body.heliocentrictrueecliptic.icrs
        else:
            # SkyCoord frame to extract lat/lon
            targ = targ.barycentrictrueecliptic
            body = body.barycentrictrueecliptic
            
        # Extract data for map plot
        #   throughout this function: we use x for longitude/ra and y for latitude/dec
        # target coords - all
        xT = getattr(targ, xattr_name).to('deg').value
        yT = getattr(targ, yattr_name).to('deg').value
        # body coords
        xB = getattr(body, xattr_name).to('deg').value
        yB = getattr(body, yattr_name).to('deg').value
        # planet names for labels
        #   now using L for Luna (moon) due to over-use of M
        names  = ['S',     'L', 'E', 'Mer', 'Ven', 'Mar', 'Jup', 'Sat', 'Ura', 'Nep', 'Plu']
        # colors = ['orange','g', 'b', 'c',   'c',   'c',   'c',   'c',   'c',   'c',   'c']
        colors = ['orange'] * len(names)
        nBodies = np.size(names)

        # time, earth-to-observatory distance, moon-to-observatory distance
        if fp_position:
            p_unit = 'au'
            fp_position.write('%s,%.6g,%.6g\n' % (time_iso,
                                       body[2].distance.to(p_unit).value,
                                       body[1].distance.to(p_unit).value))
            fp_position.close()

        # index of the "system" used by the starshade
        # (if no starshade, i.e. OS.haveOcculter == False, doesn't matter)
        # FIXME: if starshade-only, koangleArr is 1x?x?, if hybrid, is 2x?x?
        # shade_sys = 1 if (koangleArr.shape[0] == 2) else 0
        shade_sys = -1 # instead of if(...), take the last entry!
        
        if out_frames or out_movie:
            # create figure
            plt.clf()
            ax = fig.gca()

            # Note, here and elsewhere we use zorder (drawing order), largely
            # to ensure that the overlays come out on top of the patch indicating
            # the starshade exclusion zone.  So, with a starshade, the plot order
            # is: axis background, starshade patch, text/graphic overlays
            #   zorder = 30 for text and legend (drawn last)
            #   zorder = 20 +/-2 for detections + chars
            #   zorder = 10 for background graphics (unobserved stars, ecliptic line)
            #   zorder = 0 (implicit) for keepout patch

            # faint line at equator/ecliptic
            plt.hlines(0, 0, 360, colors='gray', linestyles='dashed', linewidth=0.5, zorder=10)

            for j in range(nBodies):
                # do not put any planets on the final frame, because
                # the final frame summarizes the whole mission and planet positions
                # are ephemeral
                if final_frame and (j > 0): continue

                # starshade is largely drawn as transparent, other bodies are not
                alpha = 0.17 if j <= 2 else 0.9
                # special case: if occulter, show the (large) keepout-region due to Sun (j=0)
                if OS.haveOcculter and (j == 0) and koangleArr[shade_sys,j,1] < 180:
                    # [1,j,1] => first index is for system (shade = 1), second for body (j=0 for Sun),
                    # last is for the max angle (0 for min KO, 1 for max KO)
                    koMax = koangleArr[shade_sys,j,1]
                    # lola = longitude and latitude path of circle of size koMax around xB,yB
                    lola, closed = spherical_cap(xB[j], yB[j], koMax)
                    opts = dict(facecolors='white', alpha=0.9)  # fill with white
                    #opts = dict(alpha=0.9)  # fill with white
                    paths = mplc.LineCollection(lola, linewidths=1.5, color='gray', **opts)
                    ax.add_artist(paths)

                # draw a non-circular boundary for a planet
                # lola = longitude and latitude path of circle of size koangle around xB,yB
                lola, closed = spherical_cap(xB[j], yB[j], koangleArr[0,j,0])
                opts = dict(facecolors='gray', alpha=alpha) if closed else {}
                #opts = dict(alpha=alpha) if closed else {}
                paths = mplc.LineCollection(lola, linewidths=1.5, color=colors[j], **opts)
                ax.add_artist(paths)
                    
                # marker for the text name - high zorder makes it come out on top
                ax.add_artist(plt.Text(text=names[j], x=xB[j], y=yB[j], zorder=30))
                # add a marker to indicate the center for large objects
                if koangleArr[0,j,0] > 10:
                    ax.add_artist(plt.Text(text='+', x=xB[j], y=yB[j], color='gray',
                                               horizontalalignment='center',
                                               verticalalignment='center'))

            # ====  DRM-related overlays for *observations* are within this block  ====
            # was a given star plotted?
            star_shown = np.zeros((TL.nStars,), dtype=bool)
            if OI:
                # add in observations to date: last argument is z-order
                obs_points = OI.obs_summary(time_day, not final_frame, 12)
                for obs_point in obs_points:
                    obs_star, obs_color, obs_size, obs_dict = obs_point
                    ax.add_artist(patches.Circle((OI.xpos[obs_star], OI.ypos[obs_star]),
                                                     obs_size, color=obs_color, fill=True, 
                                                     **obs_dict))
                    star_shown[obs_star] = True
                    
                # plot parameters for just below
                future_slew_mult = 0.7
                future_slew_color = [0.70, 0.9, 0.60]
                leg_props = dict(numpoints=1, labelspacing=0.2, prop={'size': 6})

                # find characterizations up to present day...
                chars = OI.char_to_date(time_day)
                # ...add in arrows for these characterizations
                for j in range(len(chars)-1):
                    r0 = chars[j]
                    r1 = chars[j+1]
                    s0 = r0['star_ind']
                    s1 = r1['star_ind']
                    # color of path from char -> char
                    pcolor_decay = 0.95 if i < Ntime-1 else 1.0 # last iteration: no decay
                    pcolor = [1 - pcolor_decay**(len(chars)-2-j)] * 3
                    # special-case arrow(j->j+1), until slew(j->j+1) has begun
                    # note, this will affect only the final arrow shown, and only between 
                    # obs[j].arrival_time and obs[j+1].arrival_time - obs[j+1].slew_time
                    #
                    # first: is slew j->j+1 a "future" slew that has not yet begun?
                    future_slew = (time_day < get_slew_begin(r1))
                    # change arrow params
                    if future_slew:
                        w_mult = future_slew_mult
                        pcolor = future_slew_color
                    else:
                        w_mult = 1.0
                    # great circle path from char -> char, as a series of line segments
                    greatc = wrap_paths((OI.xpos[s0], OI.ypos[s0]), (OI.xpos[s1], OI.ypos[s1]))
                    paths = mplc.LineCollection(greatc, linewidths=1.5*w_mult, color=pcolor)
                    ax.add_artist(paths)
                    # arrow somewhere along path
                    arrow_0, arrow_del = get_arrow(greatc, 1.5)
                    if np.linalg.norm(arrow_del) == 0:
                        continue
                    ax.arrow(arrow_0[0], arrow_0[1],
                             arrow_del[0], arrow_del[1],
                             shape='full', color=pcolor, lw=0, length_includes_head=True,
                             head_width=5*w_mult, zorder=10)

                # path-arrow legend
                # (tried with arrows rather than plt.plot, but it got hard)
                line1, = plt.plot([0,1,2], label="Spectral Char. Slew: Upcoming", lw=1.5*future_slew_mult, color=future_slew_color)
                line2, = plt.plot([0,1,2], label="Spectral Char. Slew: Recent", lw=1.5, color=[0]*3)
                line3, = plt.plot([0,1,2], label="Spectral Char. Slew: Past",   lw=1.5, color=[0.5]*3)
                # create a legend to appear at lower-right
                first_legend = plt.legend(handles=[line1, line2, line3],
                                   loc='lower right', bbox_to_anchor=(1.02, -0.245),
                                   **leg_props)
                ax = plt.gca()
                ax.add_artist(first_legend)

                # legend for DRM plot
                GS = GraphicsStyle()
                demo_dets = 4 # the count to demo the dot-size for
                l = plt.legend((
                    circle_marker(markerfacecolor=GS.det_HZ, ms=GS.count2pt2(demo_dets)),
                    circle_marker(markerfacecolor=GS.det_ok, ms=GS.count2pt2(demo_dets)),
                    circle_marker(markersize=3, markerfacecolor=GS.char_ok,   mec=GS.char_mec, mew=GS.char_lw),
                    circle_marker(markersize=3, markerfacecolor=GS.char_HZ,   mec=GS.char_mec, mew=GS.char_lw),
                    circle_marker(markersize=3, markerfacecolor=GS.char_fail, mec=GS.char_mec, mew=GS.char_lw)
                    ), (
                        f'Detection ({demo_dets}, HZ)',
                        f'Detection ({demo_dets}, non-HZ)',
                        'Spectral Char. (non-HZ)', 'Spectral Char. (HZ)', 'Failed Char. (low SNR)'),
                                   loc='lower left', bbox_to_anchor=(-0.14, -0.245),
                                   **leg_props)
                l.set_zorder(30)
        
            #------  DRM overlays end  --------
            
            # set background color of plot to gray in the case of a starshade plot
            #    if starshade: background of plot is gray, within-shade is white, other keepouts are gray
            #    no starshade: main plot is white (and no within-shade area is shown, obviously)
            if OS.haveOcculter and koangleArr[shade_sys,0,1] < 180:
                axis_color_fixed = True
                ax.patch.set_facecolor((0.92, 0.92, 0.92))
            else:
                axis_color_fixed = False
                
            # overall plot styling
            plt.minorticks_on()
            plt.xlim(0, 360);
            plt.ylim(-90, 90);
            plt.xlabel('%s [deg]' % movie_xlabel, **axlabel_style)
            plt.ylabel('%s [deg]' % movie_ylabel, **axlabel_style)
            plt.title('%s $-$ MJD %s $-$ Day #%.0f' % (time_iso, time_mjd, time_day),
                          **title_style)
            plt.tick_params(axis='both', which='major', **tick_style)

            # plot stars (visible and kept-out) with plt.scatter
            # zorder makes them appear behind the detections/characterizations,
            # to avoid clutter, we don't show the underlying target if dets/chars
            scatter_props = dict(s=5, edgecolors='none', zorder=10)
            color_ok = np.array((.6, .6, .6), ndmin=2) # [visible]  gray
            color_ko = np.array((.0, .0, .0), ndmin=2) # [kept-out] black
            plt.scatter(xT[np.logical_and( kogood[0,:,0], ~star_shown)],
                        yT[np.logical_and( kogood[0,:,0], ~star_shown)], c=color_ok, **scatter_props)
            plt.scatter(xT[np.logical_and(~kogood[0,:,0], ~star_shown)],
                        yT[np.logical_and(~kogood[0,:,0], ~star_shown)], c=color_ko, **scatter_props)
            # Need to skip transparent=True if the starshade keepout is shown,
            # otherwise the gray axis background will be turned transparent.
            # the gray axis background to transparent, which conflicts
            # with the way the occluder keepout is shown.
            savefig_opts = dict() if axis_color_fixed else dict(transparent=True) 
            # image file output
            if out_frames:
                fn = os.path.join(out_frames, 'keepout_%05d.png' % i)
                plt.savefig(fn, dpi=300, **savefig_opts)
            # movie output
            if out_movie:
                movie.grab_frame()
                if i == Ntime - 1:
                    fn = os.path.splitext(out_movie)[0] + '-final.png'
                    plt.savefig(fn, dpi=300, **savefig_opts)
    
    sys.stdout.write('\n\n')

    # need to finish the movie now, because it is linked to the current figure
    if out_movie:
        print('Movie in %s' % movie.outfile)
        movie.finish()
    if out_frames:
        print("Output image sequence in `%s'." % out_frames)

    # cumulative outputs
    if out_cume:
        # all target coords
        xa = getattr(targ, xattr_name).to('deg').value
        ya = getattr(targ, yattr_name).to('deg').value
        # cumulative good observations
        cume_good_coro = evergood_coro * (1.0 / Ntime)
        cume_good_coro_s = np.sort(cume_good_coro)
        # lengths of the maximal run for each target (in dt units)
        ok_lengths_coro = run_max(kogoods_coro) * dt
        # with occulter
        cume_good_shade = evergood_shade * (1.0 / Ntime)
        cume_good_shade_s = np.sort(cume_good_shade)
        ok_lengths_shade = run_max(kogoods_shade) * dt
        # times, for labels
        t_start = currentTimes[0]
        t_end   = currentTimes[-1]
        # turmon: 2023-11-15 -- inserted these 2 lines
        #   move from (native) MJD to ISO format b/c ISO has 'date' subformat
        t_start.format = 'iso'
        t_end.format   = 'iso'
        t_start.out_subfmt = 'date'
        t_end.out_subfmt   = 'date'

        def make_observability_cume(cume_good_sorted, title, fn):
            r'''Make and save a cumulative line plot of observability.'''
            plt.clf()
            plt.step(cume_good_sorted, np.linspace(0, 1, TL.nStars))
            plt.grid(True)
            # style it
            plt.xlabel('Portion of time observable', **axlabel_style)
            plt.ylabel('Portion of targets observed', **axlabel_style)
            plt.title('\n'.join([title, '%s $-$ %s' % (t_start.iso, t_end.iso)]),
                          **title_style)
            # 0..1 and 0..nStars y-scales
            ax1 = plt.gca()
            ax1.set_ylim((0.0, 1.0))
            #ax1.grid(True)
            par1 = ax1.twinx()
            par1.set_ylim((0, TL.nStars))
            # save it
            plt.savefig(fn, dpi=300, transparent=True)

        def make_observability_map(xa, ya, color, sz, cmap, title, label_cbar, fn):
            r'''Make and save a map-format plot of observability.'''
            plt.clf()
            scatter_props = dict(cmap=cmap, edgecolors='none')
            plt.scatter(xa, ya, c=color, s=sz, **scatter_props)
            # style it
            plt.minorticks_on()
            plt.xlim(0, 360);
            plt.ylim(-90, 90);
            plt.xlabel('%s [deg]' % movie_xlabel, **axlabel_style)
            plt.ylabel('%s [deg]' % movie_ylabel, **axlabel_style)
            plt.tick_params(axis='both', which='major', **tick_style)
            plt.title('\n'.join([title, '%s $-$ %s' % (t_start.iso, t_end.iso)]), **title_style)
            plt.colorbar(label=label_cbar)
            plt.savefig(fn, dpi=300, transparent=True)

        ### 1: cumulative plot of observability (simple line plot, not a map)
        make_observability_cume(cume_good_coro_s, 'Cumulative Target Observability: Coronagraph',
                                    os.path.join(out_cume, 'cume-obs-hist-corona.png'))
        if shadeMode:
            make_observability_cume(cume_good_shade_s, 'Cumulative Target Observability: With Occulter', 
                                    os.path.join(out_cume, 'cume-obs-hist-shade.png'))

        ### 2: map-format plot of max-continuous-observability
        make_observability_map(xa, ya, ok_lengths_coro, 18, 'jet',
                                   'Maximal Observability Map: Coronagraph', 
                                   'Longest Interval [day]', 
                                   os.path.join(out_cume, 'max-obstime-map-corona.png'))
        if shadeMode:
            make_observability_map(xa, ya, ok_lengths_shade, 18, 'jet',
                                   'Maximal Observability Map: With Occulter', 
                                   'Longest Interval [day]', 
                                   os.path.join(out_cume, 'max-obstime-map-shade.png'))

        ### 3: map-format plot of observability, plain
        cmap_obs = 'plasma' if True else 'viridis' # experimenting with colormaps
        make_observability_map(xa, ya, 100*cume_good_coro, 18, cmap_obs,
                                   'Cumulative Observability Map: Coronagraph',
                                   'Observability [%]', 
                                   os.path.join(out_cume, 'cume-obs-map-corona.png'))
        if shadeMode:
            make_observability_map(xa, ya, 100*cume_good_shade, 18, cmap_obs,
                                   'Cumulative Observability Map: With Occulter',
                                   'Observability [%]', 
                                   os.path.join(out_cume, 'cume-obs-map-shade.png'))

        ### 4: map-format plot of observability, sized by magnitude
        # sort by mag, so that bright (=large) stars plot last (=topmost)
        idx = np.argsort(TL.Vmag)
        sz_sort = np.sqrt(1e4*10**(-0.4*TL.Vmag[idx])) # sqrt() compresses vmag for size plot
        # make the plot
        make_observability_map(xa[idx], ya[idx], 100*cume_good_coro[idx], sz_sort, cmap_obs,
                                   'Cumulative Observability Map: Coronagraph (Size: Vmag)', 
                                   'Observability [%]', 
                                   os.path.join(out_cume, 'cume-obs-map-vmag-corona.png'))
        if shadeMode:
            make_observability_map(xa[idx], ya[idx], 100*cume_good_shade[idx], sz_sort, cmap_obs,
                                   'Cumulative Observability Map: With Occulter (Size: Vmag)', 
                                   'Observability [%]', 
                                   os.path.join(out_cume, 'cume-obs-map-vmag-shade.png'))
        
        # dump a CSV with visit and slew information
        OI.tour_summary(args)
        # dump a CSV with some other useful values
        fn = os.path.join(out_cume, 'target-keepout-values.csv')
        np.savetxt(fn, np.transpose(np.asarray([
            xa,
            ya,
            TL.Vmag,
            cume_good_coro,
            ok_lengths_coro])),
                       header=','.join(['lon', 'lat', 'Vmag', 'obs_frac', 'obs_maxlen']),
                       delimiter=',')

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Make EXOSIMS keepout graphics.",
                                     epilog='At least one of -m or -f or -c should be given.')
    parser.add_argument('script', metavar='SCRIPT', help='json script')
    parser.add_argument('--drm', help='file containing pickle of observations (DRM)',
                      dest='drm', metavar='FILE', default='')
    parser.add_argument('--spc', help='file of star-planet-characteristics (spc)',
                      dest='spc', metavar='FILE', default='')
    parser.add_argument('-x', '--xspecs', help='extra json spec-file to add on top of SCRIPT',
                      dest='xspecs', metavar='FILE', default='')
    parser.add_argument('-c', '--cume', help='directory for output of cumulative info',
                      dest='out_cume', metavar='DIR', default='')
    parser.add_argument('-f', '--frames', help='directory for output of still image frames',
                      dest='out_frames', metavar='DIR', default='')
    parser.add_argument('-m', '--movie', help='file for output movie',
                      dest='out_movie', metavar='FILE', default='')
    parser.add_argument('-s', '--start', help='start time [mjd], default = %(default).0f [mission start]',
                      type=float, dest='t0', default=0)
    parser.add_argument('-l', '--length', help='length [years], default = %(default).2f',
                      type=float, dest='years', default=0.5)
    parser.add_argument('-d', '--delta', help='delta-time [days], default = %(default).2f',
                      type=float, dest='delta_t', default=5.0)
    parser.add_argument('-e', '--equatorial', help='equatorial coordinates, default = False',
                      action='store_true')
    
    args = parser.parse_args()
    args.progname = os.path.basename(sys.argv[0])

    # announce how we were called, for reproducibility
    print(('%s: Invoked as: %s' % (args.progname, ' '.join(sys.argv))))

    # set umask in hopes that files/dirs will be group-writable
    os.umask(0o002)

    # plot movie in (equatorial) ra/dec coordinates vs. (ecliptic) lat/lon coordinates
    args.ra_dec_coord = not args.equatorial

    # position csv directory - not used at present
    if False:
        args.out_position = '/tmp/ko_pos'
    else:
        args.out_position = None

    if not args.out_movie and not args.out_frames and not args.out_cume:
        print('Reminder: Need one of -m or -f or -c in order to produce output.')
        # sys.exit(1)
    
    if args.ra_dec_coord:
        print('WARNING: outputs in poorly-tested ra/dec coords, use -e for equatorial.')

    # load extra specs, if any
    #    -- this isn't used (11/2023), but leaving it for now
    # argument is a filename, or if argument begins with "!", a literal.
    if args.xspecs and not args.xspecs.startswith('!'):
        assert os.path.isfile(args.xspecs), "%s is not a file." % args.xspecs
        try:
            script = open(args.xspecs).read()
            xspecs = json.loads(script)
        except ValueError as err:
            print("Error: Input xspec file `%s' improperly formatted." % (args.xspecs,))
            print("Error: JSON error was: ", err)
            # re-raise here to suppress the rest of the backtrace.
            # it is only confusing details about the bowels of json.loads()
            raise ValueError(err)
        except:
            print("Error: %s", (sys.exc_info()[0],))
            raise
    elif args.xspecs and args.xspecs.startswith('!'):
        xspecs = json.loads(args.xspecs[1:])
    else:
        xspecs = {}

    # load info about observations (DRM) from given location
    if args.drm:
        args.OI = ObserveInfo(args.drm, args.spc)
        print(args.OI.summary())
    else:
        args.OI = None
        print('Not plotting DRM.')

    make_graphics(args, xspecs)


