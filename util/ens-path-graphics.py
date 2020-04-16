#!/usr/bin/env python
#
# ens-path-graphics: Produce diagnostic summaries of sets of DRMs, including:
#   - lon/lat maps of number-of-detections
#   - lon/lat maps of number-of-slews
#   - adjacency-matrix-style heatmap plots of number of slews
#   - 3D plots (rather than lon/lat plots) of number-of-slews, illustrated with movies
#   - csv-file summaries of stars visited
# Note that this operates on sets of DRMs, not just one DRM.
#
# Usage:
#   ens-path-graphics [-a] [--outfile FILE] SCRIPT DRMs
#
# where:
#   SCRIPT is the JSON script used to generate the DRMs
#   DRMs is a list of DRM pickle files
# and, optionally:
#   -a -> skip the animation of the ensemble path 3D plot
#   --outfile FILE -> outputs filenames generated from a template,
#         which should contain 2 appearances of %s (e.g., dir/ens-path-%s.%s)

# plots summarizing DRM-sets
# turmon
#   dec 2017 (as plot-exocat)
#   mar 2018
#   sep 2018, light updates to prevent exceptions when there were 0 slews

from __future__ import division
from __future__ import print_function
import argparse
import sys
import glob
import gc
import os
import math
import six.moves.cPickle as pickle
import numpy as np
import astropy.units as u
import matplotlib as mpl; mpl.use('Agg') # not interactive: don't use X backend
import matplotlib.pyplot as plt
from matplotlib import collections as mplc
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D

import EXOSIMS.StarCatalog.EXOCAT1 as EXOCAT1
from six.moves import range
from six.moves import zip

##
## Globals
##

# filter to characterization-only?
# should perhaps be an argument
MODE = 'char'
#MODE = 'det'

# unpickling python2/numpy pickles within python3 requires this
PICKLE_ARGS = {} if sys.version_info.major < 3 else {'encoding': 'latin1'}

# use great circle paths on path plot?
GREAT_CIRCLE_PATHS = True

# better layout when there are big x-axis labels
mpl.rcParams.update({'figure.autolayout': True})

# graphics styles
axlabel_style = dict(fontsize=10, fontweight='bold')
title_style = dict(fontsize=9, fontweight='bold')
tick_style = dict(labelsize=8)
tick_tight_style = dict(labelsize=5) # tight ticks for adj matrix

# graphics constants
#cmap0 = plt.cm.get_cmap(name='viridis', lut=1024)
cmap0 = plt.cm.get_cmap(name='plasma', lut=1024) # yellow -> red -> violet
cmap0.colors = np.flip(cmap0.colors, 0)
cmap0.colors[0,0:3] = 0.7 # bottom color is gray
color_contrast = [0.2,0.8,0.2] # contrasting color for lines, etc. (green)


##
## name-related
##
# Hipparcos name, as in EXOCAT, to more common display name
HIPNAME_to_DISPLAYNAME = {
    'HIP 8102': 'Tau Ceti',
    'HIP 108870': 'Eps Indi',
    'HIP 16537': 'Eps Eridani',
    'HIP 37279': 'Procyon',
    'HIP 104214': '61 Cygni A', # 61 Cyg A
    'HIP 104217': '61 Cygni B ', # 61 Cyg B: trailing space => don't display in plots
    }
# common stars to always display in maps, if present
DISPNAME_show = ['Tau Ceti', 'Eps Indi', 'Eps Eridani', 'Procyon', '61 Cygni A']

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

def name_to_displayname(in_names, translation):
    r'''Translate each string in in_names through a table in translation.'''
    out_names = np.empty_like(in_names)
    for i in range(len(in_names)):
        out_names[i] = translation.get(in_names[i], in_names[i])
    return out_names

def name_to_index(name_list, all_names):
    r'''Return a list of indexes corresponding to appearances of name_list in all_names.'''
    out_index = []
    for name in name_list:
        # expect 1 or 0 matches per iteration here
        out_index.extend(np.where(all_names == name)[0].tolist())
    return out_index

##
## path-related
##
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
    r"""GCP = great circle points from (lon1,lat1) -> (lon2,lat2)."""
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


# container class for DRM summary info
class DRMSummary(object):
    def __init__(self, filepat, args):
        '''Load a bunch of DRMs.'''
        # 1: Load the DRMs
        drms = []
        spc = None
        Nstar = 0
        # special-case quoted wildcards for now
        if len(filepat) == 1 and '*' in filepat[0]:
            filepat = glob.glob(filepat)
        # loop over all passed-in files
        for f in filepat:
            if f.endswith('.spc'): continue
            # disabling gc during object construction speeds up by ~30% (12/2017, py 2.7.14)
            gc.disable()
            drm = pickle.load(open(f, 'rb'), **PICKLE_ARGS)
            gc.enable()
            # allow to filter the events
            if 'char' in MODE:
                drm_filter = [d for d in drm if 'char_status' in d]
                drm = drm_filter
            elif 'det' in MODE:
                drm_filter = [d for d in drm if 'det_status' in d]
                drm = drm_filter
            drms.append(drm)
            # load a spc file at first chance
            if spc is None:
                g = f.replace('pkl', 'spc').replace('/drm/', '/spc/')
                if os.path.isfile(g):
                    spc = pickle.load(open(g, 'rb'), **PICKLE_ARGS)
                    Nstar = len(spc['Name'])
        if spc is None and len(drms) > 0:
            raise ValueError('Could not find a .spc file to match the DRMs')
        if 'Name' in spc:
            spc['Name'] = np_force_string(spc['Name'])
        # save the DRMs
        self.drms = drms
        self.spc = spc
        self.Nstar = Nstar

        # 2: Find star-to-position correspondence
        # instantiate an EXOCAT object so we can harvest its data
        SC = EXOCAT1.EXOCAT1()
        # find the EXOCAT1 star corresponding to each spc star
        sc_index = np.zeros((Nstar,), dtype=np.int)
        for i in range(Nstar):
            sc_index[i] = np.where(SC.Name == spc['Name'][i])[0]
        # extract just the coordinates we care about...
        # plus, convert to heliocentric lat/lon
        self.spc_coord = SC.coords[sc_index].heliocentrictrueecliptic

    def summarize(self, args):
        '''Compute summary statistics.'''
        Nstar = self.Nstar
        Ndrm = len(self.drms)
        # list of all stars visited
        self.visited = np.array([d['star_ind'] for drm in self.drms for d in drm], dtype=np.int)
        # total cumulative number-of-visits, by star
        self.vcount = np.bincount(self.visited, minlength=Nstar)

        # slew adjacency matrix
        #   this is probably inefficient, but it may not matter
        self.visit2 = np.zeros((Nstar, Nstar), dtype=np.int)
        for drm in self.drms:
            # if only one thing in DRM, there was no slew
            if len(drm) < 2: continue
            for d1, d2 in zip(drm[0:-1], drm[1:]):
                self.visit2[d1['star_ind'], d2['star_ind']] += 1
    
        # import pdb; pdb.set_trace()
        
        # dump some info - should use a csv-writer
        if args.outfile:
            # visit counts
            fn = args.outfile % ('visits', 'csv')
            with open(fn, 'w') as f:
                arr = ['name', 'visit_mean', 'lon', 'lat', 'dist']
                f.write(','.join(arr) + '\n')
                for i, ct in enumerate(self.vcount):
                    visit_mean = ct / (1.0 * Ndrm)
                    arr = [
                        '%s' % self.spc['Name'][i],
                        '%.3f' % visit_mean,
                        '%.2f' % self.spc_coord[i].lon.value,  # [deg]
                        '%.2f' % self.spc_coord[i].lat.value,  # [deg]
                        '%.2f' % self.spc_coord[i].distance.value, # [pc]
                        ]
                    f.write(','.join(arr) + '\n')
                print("Visits written to `%s'" % fn)
            # slew counts
            fn = args.outfile % ('slews', 'csv')
            with open(fn, 'w') as f:
                pos1, pos2 = np.where(self.visit2 > 0)
                arr = ['source', 'dest', 'slews']
                f.write(','.join(arr) + '\n')
                for p1, p2 in zip(pos1, pos2):
                    arr = [
                        '%d' % p1,
                        '%d' % p2,
                        '%d' % self.visit2[p1,p2],
                        ]
                    f.write(','.join(arr) + '\n')
                print("Slews written to `%s'" % fn)
                               

def make_graphics(args, drm_info):
    r'''Bundle of statements that make the graphics - not very elegant.'''

    # convenience variables
    spc_coord = drm_info.spc_coord
    vcount = drm_info.vcount
    visit2 = drm_info.visit2
    star_names = name_to_displayname(drm_info.spc['Name'], HIPNAME_to_DISPLAYNAME)

    # ensure a blank slate
    fig = plt.figure(0)
    fig.clf()
    ax = fig.add_subplot(111)
    plt.gcf().set_size_inches(7.0,3.4)

    ###
    ### 1: Map-format plot of visits to the targets
    ###

    ## 1A: plot lines between stars that were visited in series
    # a matplotlib LineCollection might be more efficient?
    visit2_top = max(visit2.ravel())
    # we plot only the topmost N slews (N=150 at present)
    try:
        visit_floor = np.sort(visit2[visit2>0])[-150] # 150'th from top
    except IndexError:
        visit_floor = 0 # there was not a 150'th element
    if visit2_top > 0 and (len(visit2.ravel()) < 200 or True):
        for i,j in zip(*np.nonzero(visit2)):
            pair = np.array((i, j))
            if visit2[i,j] <= visit_floor: continue
            # linewidth proportional to paths from i->j
            linewidth = 5 * (visit2[i,j] / visit2_top)
            if linewidth < 0.5: continue
            if GREAT_CIRCLE_PATHS:
                # great circle path from char -> char, as a series of line segments
                greatc = wrap_paths((spc_coord.lon.value[i], spc_coord.lat.value[i]),
                                    (spc_coord.lon.value[j], spc_coord.lat.value[j]))
                paths = mplc.LineCollection(greatc, linewidths=linewidth, alpha=0.5,
                                            color=color_contrast, zorder=0)
                ax.add_artist(paths)
            else:
                ax.plot(spc_coord.lon.value[pair], spc_coord.lat.value[pair], color=color_contrast,
                        solid_capstyle='round',
                        alpha=0.5, linewidth=linewidth, zorder=0)

    # 1B: overlay a scatter plot of targets visited
    #  make the points plot in order of visit-count, so that highly-visited stars
    #  are not covered by unvisited stars
    inx = np.argsort(vcount)
    # indexes to appear with text labels
    to_label = set(name_to_index(DISPNAME_show, star_names)).union(set(inx[-6:]))

    # marker size on scatter plot: smaller if no visits
    sz = np.zeros(spc_coord.lon.value.shape) + 10.0
    sz[vcount == 0] = 2.0

    # map-format scatter plot, one point per star
    # the max(vcount, 1) prevents singularities when #visits = 0
    s = ax.scatter(spc_coord.lon.value[inx], spc_coord.lat.value[inx], c=vcount[inx], s=sz[inx],
                    norm=mpl.colors.PowerNorm(0.5, vmin=0, vmax=max(vcount.max(), 1)),
                    cmap=cmap0, zorder=10)

    # add text labels to the most-visited stars
    for i in to_label:
        if star_names[i].endswith(' '): continue # code to skip it: too close
        ax.text(spc_coord.lon.value[i], spc_coord.lat.value[i], star_names[i],
                    fontsize=4, zorder=20)

    # 1end: Annotations
    cb = plt.colorbar(s)
    cb.set_label('Expected Number of Obs. [count]')
    plt.title('Cumulative Star Paths and Visits Over %d DRMs\n%s' % (len(drm_info.drms), args.exp_name),
                  **title_style)
    plt.xlabel('Longitude [deg]', **axlabel_style)
    plt.ylabel('Latitude [deg]', **axlabel_style)
    plt.minorticks_on()
    plt.tick_params(axis='both', which='major', **tick_style)
    plt.xlim(0, 360);
    plt.ylim(-90, 90);

    # make aspect ratio close to 2 wide x 1 tall for a lat/lon plot
    plt.gcf().set_size_inches(8,3.8)
    plt.show(block=False)
    if args.outfile:
        plt.savefig(args.outfile % ('map', 'png'), dpi=300, transparent=True)

    ###
    ### 2: 3D plot of visits, with output as an animation
    ###
    scaled = True
    Ndrm = len(drm_info.drms)
    # lat, lon for convenience
    lon1 = np.deg2rad(spc_coord.lon.value)
    lat1 = np.deg2rad(spc_coord.lat.value)

    # place stars on a sphere
    if scaled:
        r = spc_coord.distance.to(u.pc).value
    else:
        r = np.ones(lon1.shape)
    # elementwise multiplication
    x1 = r * np.cos(lat1) * np.cos(lon1)
    y1 = r * np.cos(lat1) * np.sin(lon1)
    z1 = r * np.sin(lat1)

    # 2A: lines indicating slews
    fig = plt.figure(1)
    ax = fig.add_subplot(111, projection='3d')

    for i,j in zip(*np.nonzero(visit2)):
        pair = np.array((i, j))
        linewidth = 5 * (visit2[i,j] / visit2_top)
        # if linewidth < 0.2: continue
        if linewidth < 1: continue
        ax.plot(x1[pair], y1[pair], z1[pair], color=color_contrast,
                    solid_capstyle='round',
                    alpha=0.5, linewidth=linewidth, zorder=1)

    # line and dot along the N-S pole, for orientation
    ax.plot(np.array([0, 0]), np.array([0, 0]), np.array([-1, 1])*np.median(r),
                solid_capstyle='round',
                alpha=0.5, color='k', linewidth=2, zorder=0)
    ax.plot(np.array([0]), np.array([0]), np.array([0]), '.', 
                alpha=0.5, color='k', markersize=10, zorder=0)

    # 2B: scatter plot of targets visited
    #  make the points plot in order of visit-count, so that highly-visited stars
    #  are not covered by unvisited stars
    inx = np.argsort(vcount)
    # indexes to appear with text labels
    to_label = set(name_to_index(DISPNAME_show, star_names)).union(set(inx[-4:]))

    # marker size on scatter plot: smaller if no visits
    sz = np.zeros(spc_coord.lon.value.shape) + 10.0
    sz[vcount == 0] = 2.0

    s = ax.scatter(x1[inx], y1[inx], z1[inx], c=(100.0*vcount[inx]/Ndrm), s=sz[inx],
                    norm=mpl.colors.PowerNorm(0.5, vmin=0, vmax=max(vcount.max(), 1.0)),
                    cmap=cmap0, zorder=10)
    cb = plt.colorbar(s, shrink=0.5)
    cb.set_label('Expected Number of Obs. [count]')

    # add text labels to the most-visited stars
    for i in to_label:
        if star_names[i].endswith(' '): continue # code to skip it: too close
        ax.text(x1[i], y1[i], z1[i], star_names[i], fontsize=4, zorder=20)

    # Unless scaled, kill the numerical tick labels
    if scaled:
        plt.xlabel('Distance [pc]', **axlabel_style)
    else:
        ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_zticklabels([])
    plt.title('Typical Star-to-Star Slews Over %d DRMs\n%s' % (Ndrm, args.exp_name),
                  **title_style)

    # 2C: Perform the animation
    fig.show(0)

    def init():
        return fig,

    def animate(i):
        if i % 30 == 0:
            sys.stdout.write('.')
            sys.stdout.flush()
        ax.view_init(elev=min(90.0, 10.0 + i/3.0), azim=(90-i))
        return fig,

    if args.animate:
        # define the animation
        #   could not get blit=True to work
        anim = animation.FuncAnimation(fig, animate, init_func=init,
                                       frames=270, interval=20, blit=False)
        # save the animation
        if args.outfile:
            anim.save(args.outfile % ('spin', 'mp4'),
                    fps=10, extra_args=['-vcodec', 'libx264'])
        print('Animation finished.')
    fig.clf()


    ###
    ### Heat-map plot of adjacency matrix of slews
    ###
    fig = plt.figure(2)
    # indexes of all stars we slewed from or to
    #  this is a small subset of the stars visited, in the case where
    #  many DRMs have length 1
    slew_stars = np.unique(np.concatenate(np.nonzero(visit2)))

    for perm_attrib, perm_title in [('lat', 'Latitude'), ('lon', 'Longitude')]:
        fig.clf()
        ax = fig.add_subplot(111)

        # indexes of all stars visited, in permuted order
        attr = getattr(spc_coord[slew_stars], perm_attrib).value
        perm = np.argsort(attr)
        slew_stars_perm = slew_stars[perm]
        Nslew = len(slew_stars_perm)

        # plot just the submatrix of slewed stars
        # max(visit2, 1) prevents a singular transformation when there are no visits
        s = plt.imshow(visit2[np.ix_(slew_stars_perm, slew_stars_perm)], cmap=cmap0,
                   norm=mpl.colors.PowerNorm(0.5, vmin=0, vmax=max(visit2.ravel().max(), 1)))
        # plot the diagonal, to aid orientation
        plt.plot(np.arange(Nslew), np.arange(Nslew), color=color_contrast)
        if Nslew < 50:
            ax.grid(color=color_contrast, linestyle='--')
        else:
            ax.grid(False)

        #import pdb; pdb.set_trace()
        
        # indicate the star names on axes
        if Nslew > 0:
            yticklabels = [r'%s [%.0f$\!^\circ\!$]' % t for t in zip(star_names[slew_stars_perm], attr[perm])]
            plt.xticks(list(range(Nslew)), star_names[slew_stars_perm], rotation='vertical')
            plt.yticks(list(range(Nslew)), yticklabels)
            plt.tick_params(axis='both', which='major', **tick_tight_style)
        plt.title('Frequent Slews: Stars in %s Order: %d DRMs\n%s' %
                      (perm_title, Ndrm, args.exp_name), **title_style)
        plt.xlabel('Destination Star', **axlabel_style)
        plt.ylabel('Initial Star', **axlabel_style)
        # colorbar
        cb = plt.colorbar(s, shrink=0.8)
        cb.set_label('Number of Slews in DRM Set [count]')

        # make aspect ratio closer to 1:1
        fig.set_size_inches(8,7)
        # sigh, increase room at bottom + left to allow for star names
        fig.subplots_adjust(bottom=0.2, left=0.2)

        fig.show(0)
        if args.outfile:
            fig.savefig(args.outfile % ('adjacency' + '-' + perm_attrib, 'png'),
                            dpi=300, transparent=True)


def main(args):
    print('%s: beginning.' % args.progname)
    print('Loading %d file patterns.' % len(args.infile))
    drmI = DRMSummary(args.infile, args)
    print('Summarizing.')
    drmI.summarize(args)
    print('Making graphics.')
    make_graphics(args, drmI)
    print('Plots in', args.outfile % ('*', '*'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Summarize DRMs.",
                                     epilog='')
    parser.add_argument('script', metavar='SCRIPT', default='',
                            help='json script file')
    parser.add_argument('drm', metavar='DRM', nargs='*', default='.',
                            help='drm file, or directory thereof')
    parser.add_argument('--outfile', type=str, default='/tmp/exosims-tour-%s.%s',
                            help='Output file template.')
    parser.add_argument('-a', default=True, action='store_false', 
                            dest='animate', help='Skip the animation.')
    parser.add_argument('-x', '--xspecs', help='extra json spec-file to add on top of SCRIPT',
                      dest='xspecs', metavar='FILE', default='')
    args = parser.parse_args()
    
    # set umask in hopes that files/dirs will be group-writable
    os.umask(0o002)

    # do it like this for now
    args.infile = args.drm
    # args.infile = 'drms/%s.pkl' % '*'

    # name of current script, for messages
    args.progname = os.path.basename(sys.argv[0])

    if not os.path.exists(args.script):
        raise IOError('Given script file (%s) does not exist.' % args.script)

    # insert template variables in outfile if not already present
    if '%s' not in args.outfile:
        args.outfile = args.outfile + '-%s.%s'
    # ensure the directory
    outdir = os.path.dirname(args.outfile % ('test', 'test'))
    if not os.path.isdir(outdir):
        os.makedirs(outdir, 0o775)

    # not presently used
    args.matrix_order_attr = 'lon'

    # experiment name, by convention, from script
    args.exp_name = os.path.splitext(os.path.basename(args.script))[0]

    main(args)
    sys.exit(0)
