#!/usr/bin/env python
"""plot-keepout-and-obs: plot time-trace of keepout and observations within a single DRM

Usage:

  plot-keepout-and-obs.py [-o OUTPATH] SCRIPT DRM

where:

  -o OUTPATH gives an output path for results, containing two %s slots  

Most helpful Sandbox usage:

  PYTHONPATH=EXOSIMS plot-keepout-and-obs.py Scripts/FOO.json sims/FOO/drm/SEED.pkl 

where FOO.json is a script, and SEED.pkl is a DRM.  Output will be placed in the working directory 
unless -o path/to/output/%s.%s or the like is given.

Note: This imports EXOSIMS and instantiates an object based on the given SCRIPT.

Michael Turmon, JPL, 04/2019 -- created based on an idea by Dean Keithly

"""

from __future__ import print_function
import os
import sys
import csv
import six.moves.cPickle as pickle
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

# alas, this needs Exosims
import EXOSIMS
from EXOSIMS.util.get_module import get_module_from_specs

########################################
###
###  Static Data
###
########################################

SAVEFIG_OPTS = dict(dpi=300)

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
        self.specs = outspec
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
        mode_det = [mode for mode in outspec['observingModes'] if 'detection' in mode.keys()]
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
        self.ohTime = self.specs.get('ohTime', 0.2) # [days] -- overhead time
        self.settlingTime = self.specs['settlingTime'] # [days] -- settling time
        self.missionLife = self.specs['missionLife'] * 365.25 # [days]
        self.charMargin = self.specs.get('charMargin', 0.15) # [dimensionless] - Exosims default

    def get_koangles(self, OS):
        r'''Get koangles array from information in the OpticalSystem.
        The (un-pythonic) code here is copied from Prototypes/SurveySimulation so that 
        we generate exactly the koangles as expected by the caching mechanism.'''
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

    def get_keepout(self):
        r'''Instantiate Exosims object to get 2d keepout matrix and place in object.'''
        ## unused -- the below would create a whole MissionSim
        ## sim = EXOSIMS.MissionSim.MissionSim(specfile, nopar=True)
        # create just the Exosims objects we need
        specs = self.specs
        obs = get_module_from_specs(specs, 'Observatory')(**specs)
        TL  = get_module_from_specs(specs, 'TargetList')( **specs)
        TK  = get_module_from_specs(specs, 'TimeKeeping')(**specs)

        # get keepout map for entire mission
        startTime = TK.missionStart.copy()
        endTime   = TK.missionFinishAbs.copy()
        koangles = self.get_koangles(TL.OpticalSystem)
        # this will typically be cached
        #   koMap: [Mode][Target][Time], e.g. shape = (2, 405, 1827)
        #   time is in 1-day increments
        #   True = observable; False = obstructed
        koMap, koTimes = obs.generate_koMap(TL, startTime, endTime, koangles)

        # save in sim object
        self.startTime = startTime
        self.endTime = endTime
        self.koMap = koMap
        self.koTimesObj = koTimes
        # koTimes is an offset from mission start, in days
        self.koTimes = (koTimes - startTime).value
        
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
        det_t0,  det_dt,  det_sind  = [], [], []
        char_t0, char_dt, char_sind = [], [], []
        slew_t0, slew_dt, slew_sind = [], [], []
        for obs in self.drm:
            arrival_time = strip_units(obs['arrival_time'])
            star_ind = obs['star_ind']
            if ('det_info' in obs) or ('det_time' in obs):
                # a detection
                obs_time = strip_units(obs['det_time'])*timeMultiplier + ohTime + settlingTime
                det_t0.append(arrival_time)
                det_dt.append(obs_time)
                det_sind.append(star_ind)
            if has_char_info(obs):
                # a characterization -- formerly, just an "else" -- assuming that all obs
                #   were detections or characterizations
                # Note: charMargin is already included in char_time
                # TODO: check if timeMultiplier is included in char_time (not critical, we don't use it)
                # TODO: it can be that char_time = 0, perhaps not include such obs at all?
                char_t0.append(arrival_time)
                char_dt.append(strip_units(get_char_time(obs)))
                char_sind.append(star_ind)
            # 2019-09-09: relocated outside the if (char) clause for starshade-only missions
            #   which have both char and det info for some observations, and which have slews
            #   in advance of detections as well as chars.  This should handle all cases.
            # slew info -- coronagraph-only chars do not have a slew_time key
            slew_time = strip_units(obs.get('slew_time', 0.0)) # [days]
            if slew_time > 0:
                # e.g., first slew is length 0
                slew_t0.append(arrival_time-slew_time)
                slew_dt.append(slew_time)
                slew_sind.append(star_ind)
        
        self.det_t0    = np.array(det_t0)
        self.det_dt    = np.array(det_dt)
        self.det_sind  = np.array(det_sind)
        self.char_t0   = np.array(char_t0)
        self.char_dt   = np.array(char_dt)
        self.char_sind = np.array(char_sind)
        self.slew_t0   = np.array(slew_t0)
        self.slew_dt   = np.array(slew_dt)
        self.slew_sind = np.array(slew_sind)


########################################
###
###  Plot Container Class
###
########################################

class plotKeepoutContainer(object):
    """A container class for methods relating to plotting keepout and DRMs
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
        r'''Deduce a legend from what we plot.'''
        texts = []
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
        return '\n'.join(texts)

    def plot_keepout_panel(self, sim, char_only=False):
        r'''Plot one panel; char_only means only show characterization events.'''
        # plot attributes for a "collection" of events
        LF = 2 if char_only else 1
        plot_attributes = [
            ('Coronagraph', dict(linewidth=1*LF, color='blue'),  ('det_t0',  'det_dt',  'det_sind' )),
            ('Starshade',   dict(linewidth=2*LF, color='green'), ('char_t0', 'char_dt', 'char_sind' ))]
            
        # sind_show = the segment of the target list to show
        N_sind_all = sim.koMap.shape[1]
        sind_show = np.zeros(N_sind_all, dtype=bool) # False - show nothing
        if not char_only:
            sind_show[sim.det_sind] = True
        sind_show[sim.char_sind] = True
        # layer of koMap to show -- shape[0] = 1 if coronagraph-only, or = 2 for shade.
        # thus: result = 0 if detection-mode, or if coronagraph-only;
        # result = 1 if char_only and starshade
        koMap_layer = (sim.koMap.shape[0]-1) if char_only else 0
        # one figure for whole plot
        fig = plt.figure(figsize=(18,20))
        # subplots starting at t_first days, each covering t_cover days
        t_first = 0 # [days]
        t_cover = 365 # [days]
        for n_plot in range(5):
            ax = fig.add_subplot(5, 1, n_plot+1)
            t_str = 'Year %d' % (n_plot+1)
            self.plot_keepout_range(ax, plot_attributes, sim, t_str, koMap_layer, sind_show, 
                                    (t_first+n_plot*t_cover, t_first+(n_plot+1)*t_cover))

        plt.subplots_adjust(hspace=0.5)
        # add text
        if False:
            text = self.panel_legend(plot_attributes, sim)
            plt.figtext(0.02, 0.05, text, verticalalignment='center', weight='bold', fontsize=16)
        # overall title
        title = (
            'Timeline: Keepout Overlaid with Observations',
            'Showing Characterization Targets and Keepout' if char_only else \
              'Showing All Observations with Detection Keepout',
            'Detections: Blue; Characterizations: Green; Obscured: Gray')
        plt.suptitle('\n'.join(title), weight='bold', fontsize=20)
        # show the plot
        plt.show(block=False)
        fname = 'obs-keepout-char' if char_only else 'obs-keepout-all'
        plt.savefig(self.args.outpath % (fname, 'png'), **SAVEFIG_OPTS)
        # plt.savefig(self.args.outpath % (fname, 'pdf'))
        plt.close()


    def plot_keepout_range(self, ax, plot_attrs, sim, lbl, koMap_layer, sind_show, t_range):
        # plot parameters
        #   title, y-pos, y-width, color(s), (t0_attr, dt_attr)
        # put up the plots

        #import pdb; pdb.set_trace()

        # [1] -- the base keepout image
        ko_t_index = np.logical_and(sim.koTimes >= t_range[0], sim.koTimes <= t_range[1])
        # extent is given as a tuple indicating: (left, right, bottom, top)
        # in "lower" origin, the default extents are: (-0.5, numcols-0.5, -0.5, numrows-0.5)
        # we need to offset the time range by the beginning time, t_range[0]
        koMap_show = sim.koMap[np.ix_(np.array([koMap_layer]), sind_show, ko_t_index)][0]
        numrows, numcols = koMap_show.shape
        extent = (-0.5+t_range[0], numcols-0.5+t_range[0], -0.5, numrows-0.5)
        # colormap
        cmap = mpl.colors.ListedColormap(['lightgray','white'])
        # keepout map serves as background for observations
        # reminder: true = observable, false = kept out
        ax.imshow(koMap_show, aspect='auto', origin='lower', extent=extent, cmap=cmap)

        # [2] -- overlay observations
        for plot_num in range(len(plot_attrs)):
            # unpack
            name, pl_kwargs, _ = plot_attrs[plot_num]
            t0name, dtname, sindname = _
            # select obs within trange
            t0_all = getattr(sim, t0name)
            dt_all = getattr(sim, dtname)
            sind_all = getattr(sim, sindname)
            t_index = np.logical_and((t0_all+dt_all) >= t_range[0], t0_all <= t_range[1])
            self.plot_keepout_single(ax, pl_kwargs, sind_show, 
                                         t0_all[t_index], dt_all[t_index], sind_all[t_index])

        ## plot styling
        self.style_plot()
        # plot labeling info
        #y_lbl = [attrs[0] for attrs in plot_attrs if not attrs[0].startswith('-')]
        #y_pos = [attrs[1] for attrs in plot_attrs if not attrs[0].startswith('-')]
        #y_del = max([attrs[2] for attrs in plot_attrs if not attrs[0].startswith('-')])/2.0 # half-height
        ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))
        # plot labeling
        ax.set_xlim(t_range)
        #ax.set_ylim((min(y_pos)-y_del-1, max(y_pos)+y_del+1)) # pad by 1
        ax.xaxis.set_tick_params(labelsize=16)
        #ax.set_yticks(y_pos)
        #ax.set_yticklabels(y_lbl, fontsize=16)
        ax.set_ylabel('Target Index', weight='bold', fontsize=14)
        ax.set_xlabel('Time Since Mission Start [day]', weight='bold', fontsize=16)
        if sim.show_seed:
            title = 'Mission Observation and Keepout Timeline for %s: %s' % (sim.name, lbl)
        else:
            title = 'Mission Observation and Keepout Timeline: %s' % (lbl, )
        plt.title(title, weight='bold',fontsize=18)

    def plot_keepout_single(self, ax, pl_kwargs, sind_show, t0, dt, sind):
        r'''Helper function to plot one set of bars for a single instrument.'''
        # make sind2yvalue, an array mapping sind's -> as-plotted y values
        sind_used = np.nonzero(sind_show)[0]
        sind2yvalue = np.zeros(sind_show.shape, dtype=int) - 1 # will be illegal values where not set
        sind2yvalue[sind_used] = np.arange(len(sind_used))
        for i in range(len(t0)):
            # the y value for the i'th star index
            yval = sind2yvalue[sind[i]]
            # yval < 0 happens when we plot detection-only stars, in show-char-only mode
            if yval >= 0:
                # the plot lines seem to have a border that I can't turn off
                ax.plot((t0[i], t0[i]+dt[i]), (yval, yval), **pl_kwargs)

def main(args):
    # load the simulation run
    sim = SimulationRun(args.drm, args.script)
    # reduce DRM to get observation times
    sim.get_obs_times()
    # instantiate an Exosims object to get keepout for this sim
    print('%s: Instantiating an Exosims object.' % args.progname)
    sim.get_keepout()
    print('%s: Finished with the Exosims object.' % args.progname)
    # open a plot container object
    plotter = plotKeepoutContainer(args)
    # make plots
    plotter.plot_keepout_panel(sim)
    plotter.plot_keepout_panel(sim, char_only=True)
    print('%s: Done.' % args.progname)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Plot timeline of keepout and observations reported in DRMs.", epilog='')
    parser.add_argument('script', metavar='SCRIPT', default='', help='json script file')
    parser.add_argument('drm', metavar='DRM', default='', help='drm file')
    parser.add_argument('-o', default='./%s.%s', type=str, dest='outpath', help='Output file pattern.')
    #parser.add_argument('-v', default=False, action='store_true', 
    #                        dest='verbose', help='Verbosity.')
    #parser.add_argument('-s', default=True, action='store_false', 
    #                        dest='synoptic', help='Omit synoptic plot.')
    #parser.add_argument('-d', default=True, action='store_false', 
    #                        dest='dayinthelife', help='Omit day-in-the-life plot series.')
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

    # do it
    main(args)
    sys.exit(0)

