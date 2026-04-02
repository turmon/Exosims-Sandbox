#!/usr/bin/env python
#
# rad-sma-rectangle-bin-plot.py: Radius-SMA ("Kopparapu") histogram plots
# 
# This function makes plots of quantities like yield and occurrence counts, 
# binned by radius and semi-major axis (SMA). Our colloquial name for these
# plots is "Kopparapu bin plots".
#
# Usage:
#   `rad-sma-rectangle-bin-plot.py OPTS CSV`
#
# where the argument CSV is either:
#
# +  (a) the special string 'self', allowing access to baked-in
#    tables, named by the following fields:
#
#       - sdet -- sdet rates from Dulz 2019
#       - sag13 -- sag13 rates
#       - kopparapu -- kopparapu 2018 rates
# 
# +  (b) a template for (currently 2) CSV files, of the form
#       `path/to/data/reduce-%s.csv`
#   produced by the reduction script.  We need the
#   5x3 bins from files like `reduce-radlum.csv`, and
#   earth chars from `reduce-earth.csv`.  Several specific
#   columns ("fields") can be extracted (from each).
#
# In addition:
# 
# +  -f FIELD: comma-separated name(s) of the field(s) within
#     the CSV to plot - use any or all of:
# 
#        -  char_{full,strict}
#        -  char_tput_{full,strict}
#        -  population
#     or, for the "self" target, one or more of the fields
#     noted above.
#     By default, all fields are given.
# 
# +  -t, --title : set graph title.  By default, a string
#     derived from the field name is used, so this option
#     is not much needed.
# +  --eta : suppress \eta = ... text label in table.  Generally 
# +  --sigma : add +/- sigma suffix to the numbers in the table
# +  --quantile : add ^{U}_{L} where U and L are the distance
#         to the upper and lower quantile (U = q75 - mean, etc.)
# +  -o : specify an output-file template, with *two* %s placeholders
#        for a graph type (derived from the field name)
#        and a file type (pdf/png).  By default: ./det-rad-sma-%s.%s
# +  -h : help
#
# Sample usage:
#```shell
#   rad-sma-rectangle-bin-plot.py --sigma --eta -t 'SDET Rates' -f sdet self
#   rad-sma-rectangle-bin-plot.py --sigma --eta -t 'CSV Rates' -f char_full csv/reduce-%s.csv
#```
#
# For more on usage, use the -h option.
# Some options may be described there but not documented here.

# author:
#  Michael Turmon, JPL, 2019

import os
import sys
import csv
import copy
import math
import argparse
from collections import defaultdict
from pathlib import Path
import numpy as np
# this import must work: fail fast if it doesn't
from reduce_drm_tools.PlanetBins import RpLBins
from reduce_drm_tools import utils

import matplotlib as mpl; mpl.use('Agg') # not interactive: don't use X backend
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle, Polygon
from matplotlib.ticker import FuncFormatter

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
    # r'' to shield \e and \o from interpretation
    if eta:
        symbol = r'\eta%s = ' % (r'_\oplus' if earth else '')
    else:
        symbol = ''
    txt = '%s%#.3g' % (symbol, x)
    if ranges:
        # we are allowed to show a range
        if x_lo != None and x_hi != None:
            # low/hi range
            txt += '^{%+.2g}_{%+.2g}' % (x_hi - x, x_lo - x)
        elif x_lo != None and x_hi == None:
            # standard deviation (r'' because \p would look like \n)
            # superscript \pm to make the range in smaller font
            txt += r'^{{\pm}%#.2g}' % (x_lo, )
    return '$' + txt + '$'

def make_koppa_boxes(args, ax, hist):
    r'''Render the table.'''
    ## Set up a bunch of data
    SHOW_RANGES = True
    Text_style_base = dict(fontsize=14, 
                           horizontalalignment='center',
                           verticalalignment='center')
    # eta: used within the plot area for all rates or counts
    Text_style_eta = copy.copy(Text_style_base);
    Text_bubble = dict(bbox=dict(boxstyle='round', facecolor='white', linewidth=0, alpha=0.4))

    # instantiate a binner to get planet categories
    binner = RpLBins()
    Rp_bins = binner.Rp_bins

    #Earth_style = dict(facecolor=None, edgecolor='lightgreen', linewidth=2.0, hatch='/')
    Earth_style = dict(fill=False,
                       facecolor=None,
                       edgecolor='lightgreen',
                       linewidth=2.0,
                       alpha=0.7, # tiny bit of transparency
                       hatch='/')

    # L bins are in high...low order; a-bins are in low-high order
    #   colors below tweaked repeatedly, formerly "tomato"
    style_for_L = [dict(facecolor=fc) for fc in ('xkcd:pastel red', 'dodgerblue', 'lightskyblue')]
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
            # ad hoc label bumps to make room for earth
            if r_inx == 1 and L_inx == 1:
                r_mid = r_mid * 1.12 # upward
            elif r_inx == 0 and L_inx == 1:
                r_mid = r_mid * 0.85 # downward
            info = hist[L_inx + r_inx*Lbin_num]
            txt = make_text_slug(info[0], info[1], info[2], SHOW_RANGES, eta=args.eta)
            tbox = ax.text(a_mid, r_mid, txt, **Text_style_eta)
    # Nevada-shaped polygon for Earths
    #   (zorder tweak needed for visibility)
    #   TBD: in is_earthlike(), Rp_hi (upper) is handled differently from Rp_lo (lower)
    #   If this plot is to reflect the is_earthlike() boundary, it must duplicate the
    #   is_earthlike() logic
    earth_x = np.array([binner.Earth_SMA_lo, binner.Earth_SMA_hi, binner.Earth_SMA_hi, binner.Earth_SMA_lo])
    earth_y = np.array([binner.Earth_Rp_hi,  binner.Earth_Rp_hi,  binner.Earth_Rp_lo2, binner.Earth_Rp_lo1])
    earth_gon = Polygon(np.vstack((earth_x, earth_y)).T, zorder=2, **Earth_style)
    ax.add_artist(earth_gon)
    # put in the Earth eta label
    if len(hist) > 15:
        info = hist[-1]
        txt = make_text_slug(info[0], info[1], info[2], SHOW_RANGES, earth=True, eta=args.eta)
        a_mid = np.exp(np.mean(np.log(earth_x)))
        r_mid = np.exp(np.mean(np.log(earth_y)))*0.95
        tbox = ax.text(a_mid, r_mid, txt, color='darkgreen', **Text_style_eta)
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
        # plot albedo box (Gray A_G = ... near left axis)
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

    # tidy up 0.09999999 on x-axis
    formatter = FuncFormatter(lambda y, _: '{:.8g}'.format(y))
    ax.xaxis.set_major_formatter(formatter)

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
        formatter = FuncFormatter(lambda y, _: '{:.8g}'.format(y))
        axis.set_major_formatter(formatter)

    # dump the plot
    plt.tight_layout()
    if False:
        plt.show() # OK for interactive, not for scripted
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
    os.umask(0o002)

    # load local reduction parameters
    # find enclosing directory
    # directory = Path(os.path.dirname(args.outfile % ('dummy', 'txt'))).parent
    directory = Path(os.path.dirname(args.csv))
    args.reduce_config = utils.load_reduce_config(directory, log_origin=args.progname)
    if args.reduce_config is None:
        print(f'{args.progname}: No local reduction config file ({utils.REDUCTION_CONFIG})')
        # so we can always assume it's a dict
        args.reduce_config = {}
    else:
        print(f'{args.progname}: Loaded reduction config: {args.reduce_config["_config_filename"]}')
    # customize the overall binner class
    fails = RpLBins.customize_parameters(args.reduce_config)
    if fails:
        print(f'{args.progname}: Warning: {len(fails)} unused attribute(s) in {args.reduce_config["_config_filename"]}')
        print(f'{args.progname}: Warning: Unused attributes: {", ".join(fails)}')
    elif args.reduce_config:
        print(f'{args.progname}: Note: Local customization successful.')
        
    args.field_list = args.fields.split(',')
    assert len(args.field_list) > 0, 'Need at least one field to be given (-f FIELDS)'

    print('%s: Producing %d radius/SMA plots' % (args.progname, len(args.field_list)))

    # do it
    main(args)
    sys.exit(0)


