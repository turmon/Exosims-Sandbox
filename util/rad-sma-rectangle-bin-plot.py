#!/usr/bin/env python
#
# Radius-SMA ("Kopparapu") bin plot of planet properties (characterizations, population)
#
# Usage:
#   rad-sma-rectangle-bin-plot.py OPTS CSV
# where the argument CSV is either:
#   (a) the special string 'self', allowing access to baked-in
#   tables, named by the following fields:
#       sdet -- sdet rates from Dulz 2019
#       sag13 -- sag13 rates
#       kopparapu -- kopparapu 2018 rates
#   (b) a template for (currently 2) CSV files, of the form
#       path/to/data/reduce-%s.csv
#   produced by the reduction script.  We need the
#   5x3 bins from files like "reduce-radlum.csv", and
#   earth chars from "reduce-earth.csv".  Several specific
#   columns ("fields") can be extracted (from each).
# where OPTS is:
#   -f FIELD: comma-separated name(s) of the field(s) within
#     the CSV to plot - use any or all of:
#          char_{full,strict}
#          char_tput_{full,strict}
#          population
#     or, for the "self" target, one or more of the fields
#     noted above.
#     By default, all fields are given.
#   -t, --title : set graph title.  By default, a string
#     derived from the field name is used, so this option
#     is not much needed.
#   --eta : suppress \eta = ... text label in table.  Generally 
#   --sigma : add +/- sigma suffix to the numbers in the table
#   --quantile : add ^{U}_{L} where U and L are the distance
#         to the upper and lower quantile (U = q75 - mean, etc.)
#   -o : specify an output-file template, with *two* %s placeholders
#        for a graph type (derived from the field name)
#        and a file type (pdf/png).  By default: ./det-rad-sma-%s.%s
#   -h : help
#
# Sample usage:
#   rad-sma-rectangle-bin-plot.py --sigma --eta -t 'SDET Rates' -f sdet self
#   rad-sma-rectangle-bin-plot.py --sigma --eta -t 'CSV Rates' -f char_full csv/reduce-%s.csv
#
# For more on usage, use the -h option.
# Some options may be described there but not documented here.

# author:
#  Michael Turmon, JPL, 2019

from __future__ import print_function
import os
import sys
import csv
import copy
import math
import argparse
from collections import defaultdict
import numpy as np

import matplotlib as mpl; mpl.use('Agg') # not interactive: don't use X backend
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle, Polygon

## mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
## #rc('font',**{'family':'serif','serif':['Palatino']})
## mpl.rc('text', usetex=True)

########################################
###
###  Static Data
###
########################################

SAVEFIG_OPTS = dict(dpi=300)

# Copied and wrangled from Table 3 of: Ravi Kumar Kopparapu, et al.,
# Exoplanet Classification and Yield Estimates for Direct Imaging Missions,
# ApJ, 856:122, 2018.
# (Not presently used.)
KoppaData = [
    ["Hot_Rocky",         "(182-1.0)",      "0.5-1.0",  0.22,  0.67,  2.04], 
    ["Warm_Rocky",        "(1.0-0.28)",     "0.5-1.0",  0.09,  0.30,  1.04], 
    ["Cold_Rocky",        "(0.28-0.0035)",  "0.5-1.0",  0.50,  1.92,  7.61], 
    ["Hot_Super-Earths",  "(187-1.12)",     "1.0-1.75", 0.21,  0.47,  1.04], 
    ["Warm_Super-Earths", "(1.12-0.30)",    "1.0-1.75", 0.087, 0.21,  0.54], 
    ["Cold_Super-Earths", "(0.30-0.0030)",  "1.0-1.75", 0.50,  1.42,  4.14], 
    ["Hot_Sub-Neptunes",  "(188-1.15)",     "1.75-3.5", 0.29,  0.48,  0.79], 
    ["Warm_Sub-Neptunes", "(1.15-0.32)",    "1.75-3.5", 0.12,  0.22,  0.41], 
    ["Cold_Sub-Neptunes", "(0.32-0.0.0030)","1.75-3.5", 0.77,  1.63,  3.52], 
    ["Hot_Sub-Jovians",   "(220-1.65)",     "3.5-6.0",  0.05,  0.07,  0.12], 
    ["Warm_Sub-Jovians",  "(1.65-0.45)",    "3.5-6.0",  0.04,  0.07,  0.13], 
    ["Cold_Sub-Jovians",  "(0.45-0.0030)",  "3.5-6.0",  0.58,  1.35,  3.19], 
    ["Hot_Jovians",       "(220-1.65)",     "6.0-14.3", 0.028, 0.056, 0.11], 
    ["Warm_Jovians",      "(1.65-0.40)",    "6.0-14.3", 0.023, 0.053, 0.12], 
    ["Cold_Jovians",      "(0.40-0.0025)",  "6.0-14.3", 0.34,  1.01,  3.07], 
    ["Eta Earth", 		  "",               "",      	None,  0.24,  None], # not from the table
    ]

# SAG13 data taken from report:
# * Tabulation of SDET occurrence rate distributions (right cols) with
#   older SAG13 occurrence rates (second column), for the
#   Kopparapu et al. 2018 planet type definitions.
# * Also, last row, $\eta_{Earth}$ values for limits of
# $0.95<a<1.67$ AU and  $0.8/\sqrt{a}< R_p <1.4 R_{\Earth}$.
# Columns:
#   Planet Type | SAG 13 | Optimistic | Nominal | Pessimistic
SAG13Data = [
    ["Hot rocky", 			0.67, 	1.82, 	0.64, 	0.22],
    ["Warm rocky", 			0.30, 	1.07, 	0.31, 	0.09],
    ["Cold rocky", 			1.92, 	3.80, 	1.89, 	0.50],
    ["Hot super-Earths", 	0.47, 	0.88, 	0.43, 	0.21],
    ["Warm super-Earths", 	0.21, 	0.56, 	0.22, 	0.09],
    ["Cold super-Earths", 	1.42, 	1.36, 	1.33, 	0.51],
    ["Hot sub-Neptunes", 	0.48, 	0.66, 	0.44, 	0.28],
    ["Warm sub-Neptunes", 	0.22, 	0.41, 	0.23, 	0.12],
    ["Cold sub-Neptunes", 	1.63, 	1.19, 	1.38, 	0.78],
    ["Hot sub-Jovians", 	0.07, 	0.10, 	0.07, 	0.05],
    ["Warm sub-Jovians", 	0.07, 	0.13, 	0.07, 	0.04],
    ["Cold sub-Jovians", 	1.35, 	1.14, 	1.06, 	0.58],
    ["Hot Jovians", 		0.056, 	0.07, 	0.06, 	0.05],
    ["Warm Jovians", 		0.053, 	0.13, 	0.08, 	0.06],
    ["Cold Jovians", 		1.01, 	1.48, 	0.85, 	0.45],
    ["Eta Earth", 			0.24, 	0.71, 	0.24, 	0.09],
    ]

########################################
###
###  Classes
###
########################################

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

    # extracted from the below
    Earth_SMA_lo = 0.95
    Earth_SMA_hi = 1.67
    Earth_Rp_hi = 1.4
    Earth_Rp_lo1 = 0.8/np.sqrt(Earth_SMA_lo) # left bin boundary
    Earth_Rp_lo2 = 0.8/np.sqrt(Earth_SMA_hi) # right boundary
    
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


############################################################
##
## Program-specific code below
##
############################################################

def get_static_hist(args, locator):
    r'''Given a tag, returns data, a list of 15+1 tuples, each of (value, low, high).'''
    if locator == 'sdet':
        hist = [(x[3], x[4], x[2]) for x in SAG13Data] # Dulz 2019 numbers w/ quantiles
    elif locator == 'sag13':
        hist = [(x[1], None, None) for x in SAG13Data] # no error bars
    elif locator == 'kopparapu':
        hist = [(x[4], x[3], x[5]) for x in KoppaData] # K. 2018 numbers w/ quantiles
    else:
        assert False, 'Do not recognize %s' % locator
        hist = None
    return hist
        
def load_hist(args, field):
    r'''Returns data, a list of 15+1 triples, each of (value, low, high).

    If low, high are returned as None, it's interpreted as a raw number to go into the box
    (later, when boxes are rendered).  If just high is None, it's interpreted as a standard
    error on an error bar ("sigma").  If all are numeric, a quantile (x^{U}_{L}) will be 
    rendered.'''
    fn = args.csv
    assert fn.count('%') == 1, 'Need ONE filename template (%%s) in the file component of "%s"' % fn
    # 1: grab the 3x5 binned results
    info = []
    fn_bins = fn % 'radlum'
    try:
        with open(fn_bins, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                info.append(row)
    except IOError:
        print('%s: Fatal.  Unable to open source (histogram) CSV file %s' % (args.progname, fn_bins))
        raise
    assert len(info) == 15, 'Expected 15 rows in %s' % fn_bins
    # 2: grab the Earth results
    earth = []
    fn_earth = fn % 'earth'
    try:
        with open(fn_earth, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                earth.append(row)
    except IOError:
        print('%s: Fatal.  Unable to open source (earth) CSV file %s' % (args.progname, fn_earth))
        raise
    assert len(earth) == 1, 'Expected 1 row in %s' % fn_earth
    # 3: extract the fields and compile
    hist = []
    # a: 5x3 histogram
    field_bins = 'h_RpL_%s' % field
    f1 = field_bins+'_mean'
    assert f1 in info[0], 'Correct field (%s) not found in CSV radius/luminosity file' % f1
    for row in info:
        if args.sigma:
            hist.append((float(row[field_bins+'_mean']), float(row[field_bins+'_std']), None))
        elif args.quantile:
            hist.append((
                float(row[field_bins+'_q50']),
                float(row[field_bins+'_q25']),
                float(row[field_bins+'_q75'])))
        else:
            hist.append((float(row[field_bins+'_mean']), None, None))
    # b: earth
    field_earth = 'exoE_%s' % field
    f1 = field_earth+'_mean'
    assert f1 in earth[0], 'Correct field (%s) not found in CSV earth-char file' % f1
    row = earth[0]
    if args.sigma:
        hist.append((float(row[field_earth+'_mean']), float(row[field_earth+'_std']), None))
    elif args.quantile:
        hist.append((
            float(row[field_earth+'_q50']),
            float(row[field_earth+'_q25']),
            float(row[field_earth+'_q75'])))
    else:
        hist.append((float(row[field_earth+'_mean']), None, None))
    # 4: done
    return hist
        
def make_text_slug(x, x_lo, x_hi, ranges=False, earth=False, eta=True):
    r'''Generate the "slug" of text that is placed into the cells of the 
    tabular plot.'''
    # box preface: "FOO = x", or just "x".  Not general enough, because
    # eta is not the right FOO for some fields.
    if eta:
        symbol = '\eta%s = ' % ('_\oplus' if earth else '')
    else:
        symbol = ''
    txt = '%s%.2g' % (symbol, x)
    if ranges:
        # we are allowed to show a range
        if x_lo != None and x_hi != None:
            # low/hi range
            txt += '^{%+.2g}_{%+.2g}' % (x_hi - x, x_lo - x)
        elif x_lo != None and x_hi == None:
            # standard deviation
            txt += '\pm%.2g' % (x_lo, )
    return '$' + txt + '$'

def make_koppa_boxes(args, ax, hist):
    r'''Render the table.'''
    ## Set up a bunch of data
    SHOW_RANGES = True
    Text_style_base = dict(fontsize=15, 
                            horizontalalignment='center',
                            verticalalignment='center')
    Text_style_eta = copy.copy(Text_style_base);
    Text_bubble = dict(bbox=dict(boxstyle='round', facecolor='white', linewidth=0, alpha=0.4))
    binner = RpLBins()
    Rp_bins = binner.Rp_bins

    #Earth_style = dict(facecolor=None, edgecolor='lightgreen', linewidth=2.0, hatch='/')
    Earth_style = dict(fill=False, facecolor=None, edgecolor='lightgreen', linewidth=2.0, hatch='/')

    # L bins are in high...low order; a-bins are in low-high order
    style_for_L = [dict(facecolor=fc) for fc in 'tomato', 'dodgerblue', 'lightskyblue']
    style_for_all = dict(edgecolor='white')
    x_axis_bump = [0.8, 1.0, 1.3] # visual tweak of text label locations
    
    # Create list for all the Kopparapu boxes
    # iterate to set up 5x3 boxes
    koppa_boxes = []
    # number of Rp bins
    Rbin_num = len(Rp_bins) - 1
    for r_inx in range(Rbin_num):
        r0, r1 = Rp_bins[r_inx:r_inx+2]
        L_bins = binner.L_bins[r_inx]
        # conversion: L = 1/a^2 => a = 1/sqrt(L)
        a_bins = [1.0/math.sqrt(x) for x in L_bins]
        Lbin_num = len(L_bins) - 1
        for L_inx in range(Lbin_num):
            a0, a1 = a_bins[L_inx:L_inx+2]
            rect = Rectangle((a0, r0), a1-a0, r1-r0)
            koppa_boxes.append((rect, style_for_L[L_inx]))
            # geometric mean for log-scale
            a_mid = math.sqrt(a1 * a0) * x_axis_bump[L_inx]
            r_mid = math.sqrt(r1 * r0)
            if r_inx == 1 and L_inx == 1:
                r_mid = r_mid * 1.15 # ad hoc label bump
            elif r_inx == 0 and L_inx == 1:
                r_mid = r_mid * 0.8 # ad hoc label bump
            info = hist[L_inx + r_inx*Lbin_num]
            txt = make_text_slug(info[0], info[1], info[2], SHOW_RANGES, eta=args.eta)
            tbox = ax.text(a_mid, r_mid, txt, **Text_style_eta)
    # Nevada-shaped polygon for Earths
    earth_x = np.array([binner.Earth_SMA_lo, binner.Earth_SMA_hi, binner.Earth_SMA_hi, binner.Earth_SMA_lo])
    earth_y = np.array([binner.Earth_Rp_hi,  binner.Earth_Rp_hi,  binner.Earth_Rp_lo2, binner.Earth_Rp_lo1])
    earth_gon = Polygon(np.vstack((earth_x, earth_y)).T, **Earth_style)
    ax.add_artist(earth_gon)
    # put in the Earth eta label
    if len(hist) > 15:
        info = hist[-1]
        txt = make_text_slug(info[0], info[1], info[2], SHOW_RANGES, earth=True, eta=args.eta)
        a_mid = np.exp(np.mean(np.log(earth_x)))
        r_mid = np.exp(np.mean(np.log(earth_y)))*0.95
        tbox = ax.text(a_mid, r_mid, txt, **Text_style_eta)
    else:
        print('No Earth info found.')

    # separated out only for historical reasons
    for box, style in koppa_boxes:
        # Create patch collection with specified style
        style1 = copy.copy(style_for_all)
        style1.update(style)
        pc = PatchCollection([box], **style1)
        # Add collection to axes
        ax.add_collection(pc)

    if args.albedo:
        # plot albedo box
        albedo_Rp = 1.4 # earth radii: below, is 0.2, above, is 0.5
        albedo_box = Rectangle((0, 0), 100, albedo_Rp, color='gray', zorder=-1)
        ax.add_artist(albedo_box)
        # albedo text
        albedo_sma_anchor = 0.06 # arbitrary, visually-determined
        txt = '$A_G = 0.2$'
        this_text_style = copy.copy(Text_style_base)
        this_text_style['verticalalignment'] = 'top' # label will go below line
        this_text_style['rotation'] = 90
        this_text_style.update(Text_bubble)
        tbox = ax.text(albedo_sma_anchor, albedo_Rp/1.3, txt, **this_text_style)
        txt = '$A_G = 0.5$'
        this_text_style = copy.copy(Text_style_base)
        this_text_style['verticalalignment'] = 'bottom' # label above line
        this_text_style['rotation'] = 90
        tbox = ax.text(albedo_sma_anchor, albedo_Rp*2, txt, **this_text_style)

    # seems needed to show a plot
    ax.plot(1.0, 1.0)

    ax.set_yscale('log')
    ax.set_xscale('log')
    #ax.grid(True)
    return pc

## Static data mapping fields to titles for code below

FIELD_INFO = dict(
    # built-ins
    sdet=('SDET Planet Population', 'count', 'count'),
    sag13=('SAG13 Planet Population', 'count', 'count'),
    kopparapu=('Kopparapu (2018) Planet Population', 'count', 'count'),
    # CSV
    population=('Planet Population by Type', 'count', 'count'),
    char_full=('Number of Full Characterizations', 'count', 'count'),
    char_strict=('Number of Strict Characterizations', 'count', 'count'),
    char_tput_full=('Throughput: Full Characterizations', 'ratio', '%'),
    char_tput_strict=('Throughput: Strict Characterizations', 'ratio', '%'),
    )

def field_to_scale_factor(args, field):
    try:
        unit = FIELD_INFO[field][2] if args.percent else FIELD_INFO[field][1]
        if unit == '%':
            scale = 100
        else:
            scale = 1
        return scale
    except KeyError:
        return 1
    
def field_to_title(args, field):
    # decoder ring, from field to printed title
    try:
        name = FIELD_INFO[field][0]
        unit = FIELD_INFO[field][2] if args.percent else FIELD_INFO[field][1]
        if unit:
            title = '%s [%s]' % (name, unit)
        else:
            title = name
        return title
    except KeyError:
        return ''

def get_short_name(locator):
    r'''Return a base name for plot output file; locator validation happened already.'''
    return locator.replace('_', '-')

def retrieve_data(args, field):
    # extract specific numbers to plot
    if args.csv == 'self':
        hist = get_static_hist(args, field)
    else:
        hist = load_hist(args, field)
    scale = field_to_scale_factor(args, field)
    # apply scale factor
    def apply_scale(s, x):
        return None if x is None else s*x
    hist_scaled = []
    for row in hist:
        hist_scaled.append([apply_scale(scale, x) for x in row])
    return hist_scaled

def plot_rects(args, field):
    r'''Make the plot, and dump it to a file.'''
    # Get data
    hist = retrieve_data(args, field)
    # create (A, Rp) boxes using data in hist
    fig, ax = plt.subplots(1)
    make_koppa_boxes(args, ax, hist)

    # Decorations: titles, axis labels
    ax.set_ylabel('$R_p$ [Earth radii]', fontsize=14, labelpad=-8)
    ax.set_xlabel('Semi-Major Axis [AU]', fontsize=14)
    title = field_to_title(args, field)
    if args.title:
        title = args.title
    elif args.name:
        title = '\n'.join((title, args.name))
    if title:
        ax.set_title(title, fontsize=16)
    ax.set_facecolor('lightgray')
    # ticks as numbers not exponentials
    for axis in [ax.xaxis, ax.yaxis]:
        formatter = mpl.ticker.FuncFormatter(lambda y, _: '{:.16g}'.format(y))
        axis.set_major_formatter(formatter)

    # dump the plot
    plt.tight_layout()
    plt.show()
    short_name = get_short_name(field)
    for file_type in ('png', 'pdf'):
        fn = args.outpath % (short_name, file_type)
        plt.savefig(fn, **SAVEFIG_OPTS)
    plt.close()

def main(args):
    for field in args.field_list:
        plot_rects(args, field)

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Binned plot in radius/luminosity space of planet properties.", epilog='')
    parser.add_argument('csv', metavar='csv', type=str, help='csv file (or "self")')
    parser.add_argument('-f', dest='fields', type=str, help='comma-separated field list (must be non-empty)')
    parser.add_argument('-o', default='./det-%s.%s', type=str, dest='outpath', help='Output file pattern.')
    parser.add_argument('-n', default='', type=str, dest='name', help='Explicitly-given name for plot labels.')
    parser.add_argument('-t', default='', type=str, dest='title', help='Explicitly-given full plot title.')
    parser.add_argument('--sigma', default=False, action='store_true', 
                            dest='sigma', help='Show standard deviation.')
    parser.add_argument('--percent', default=False, action='store_true', 
                            dest='percent', help='Ratios shown in percent rather than in [0,1].')
    parser.add_argument('--quantile', default=False, action='store_true', 
                            dest='quantile', help='Show quantile range (q25/q75).')
    # the below is basically required at this point, because most things plotted aren't eta!
    parser.add_argument('--eta', default=True, action='store_false', 
                            help='Omit eta = label.')
    parser.add_argument('--albedo', default=True, action='store_false', 
                            help='Omit albedo A_G = ... annotation.')
    #parser.add_argument('-v', default=False, action='store_true', 
    #                        dest='verbose', help='Verbosity.')
    args = parser.parse_args()
    args.progname = os.path.basename(sys.argv[0])
    
    # announce how we were called, for reproducibility
    print('%s: Invoked as: %s' % (args.progname, ' '.join(sys.argv)))

    if args.outpath.count('%s') != 2:
        sys.stderr.write('Fatal.  Outpath (%s) must have two %%s patterns.' % args.outpath)
        sys.exit(1)

    # set umask in hopes that files will be group-writable
    os.umask(0002)

    args.field_list = args.fields.split(',')
    assert len(args.field_list) > 0, 'Need at least one field to be given (-f FIELDS)'

    print('%s: Producing %d radius/SMA plots' % (args.progname, len(args.field_list)))

    # do it
    main(args)
    sys.exit(0)


