#!/usr/bin/env python
"""plot-star-obs-trace: plot time-trace of keepout and observations within a single DRM

Usage:

  plot-str-obs-trace.py [-o OUTPATH] DRM

where:

  -o OUTPATH gives an output path for results, containing two %s slots  

Most helpful Sandbox usage:

  plot-star-obs-trace.py sims/FOO/drm/SEED.pkl 

SEED.pkl is a DRM.  Output will be placed in the working directory 
unless -o path/to/output/%s.%s or the like is given.

If present, the corresponding SPC will be deduced and loaded to give the
correct star names. If not, a warning is given, but it's not a failure.

Note: does not need to import EXOSIMS

Michael Turmon, JPL, 07/2023

"""

from __future__ import print_function
import os
import sys
import csv
import warnings
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
#import EXOSIMS
#from EXOSIMS.util.get_module import get_module_from_specs

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
    def __init__(self, drmfile):
        """ loads pkl and script files
        Args:
            drmfile (string) - path to drm pickle file to load
        Return:
            initialized SimulationRun object

        Instantiates:
            DRM (list) - a list of observations
        """
        try:
            with open(drmfile, 'rb') as f:
                DRM = pickle.load(f, **PICKLE_ARGS)
        except:
            sys.stderr.write('Failed to open DRM file "%s"\n' % drmfile)
            raise
        # load a spc file
        load_spc = True
        if load_spc:
            g = drmfile.replace('pkl', 'spc').replace('/drm/', '/spc/')
            if os.path.isfile(g):
                spc = pickle.load(open(g, 'rb'), **PICKLE_ARGS)
            else:
                sys.stderr.write(f'No .spc file at {g}, continuing.\n')
                spc = None

        # set up object state
        #self.seed = int(os.path.splitext(os.path.basename(f))[0])
        self.spc = spc # = None, if SPC not present

        # save this in the object state
        self.drm = DRM
        self.name = os.path.splitext(os.path.basename(drmfile))[0] 
        # expedient way to suppress the numerical seed for presentation plots
        self.show_seed = True
        try:
            if os.path.exists(os.path.join(os.path.dirname(drmfile), '..', 'EnsembleName.txt')):
                self.show_seed = False
        except:
            pass


    def get_obs_trace(self):
        r'''Extract observation information, place in self object.
        '''
        # obs_trace key:
        #   = 0      ==> no det obs yet
        #   = m > 0  ==> m det obs (successful or not)
        #   = -1     ==> one fully successful char obs
        #   = -2     ==> one partially successful char obs
        #   = -3     ==> one failed char obs
        sInd_set = set(obs['star_ind'] for obs in self.drm)
        # map that converts sInd to sobsInd
        s2so = {sInd:soInd for soInd, sInd in enumerate(sorted(sInd_set))}
        Nsobs = len(s2so)
            
        # variable initializations
        obs_trace = np.zeros((Nsobs, len(self.drm)+1))
        det_ok, det_nok, det_nul, det_iwa, det_owa, det_fail = [], [], [], [], [], []
        char_ok, char_nok, char_part, char_miss = [], [], [], []
        arrival_times = []
        det_num_obs = np.zeros((Nsobs,))
        char_num_obs = np.zeros_like(det_num_obs)
        planet_num = np.zeros_like(det_num_obs)
        char_success = np.zeros_like(det_num_obs)

        # find observation trace
        for nobs, obs in enumerate(self.drm):
            arrival_times.append(obs['arrival_time'].value)
            # we are here setting nobs+1
            # copy obs_trace forward, if it shows detections only
            for sob in range(Nsobs):
                if obs_trace[sob, nobs] > 0:
                    obs_trace[sob, nobs+1] = obs_trace[sob, nobs]
                elif obs_trace[sob, nobs] < 0:
                    obs_trace[sob, nobs+1] = obs_trace[sob, nobs]
            sInd = obs['star_ind']
            soInd = s2so[sInd]
            if ('det_info' in obs) or ('det_time' in obs):
                # a detection
                obs_trace[soInd, nobs+1] = obs_trace[soInd, nobs] + 1
                det_num_obs[soInd] += 1
                planet_num[soInd] = len(obs['det_status'])
                if any(obs['det_status'] > 0):
                    det_ok.append((soInd, nobs))
                    any_det = True
                else:
                    det_nok.append((soInd, nobs)) # failure (1 of several types)
                    any_det = False
                # itemize the 3 types of fail - let any success win out
                if len(obs['det_status']) == 0:
                    det_nul.append((soInd, nobs))
                if not any_det and any(obs['det_status'] == 0):
                    det_fail.append((soInd, nobs)) # worst of the fails
                else:
                    if not any_det and any(obs['det_status'] == -2):
                        det_owa.append((soInd, nobs))
                    if not any_det and any(obs['det_status'] == -1):
                        det_iwa.append((soInd, nobs))
            if has_char_info(obs):
                char_num_obs[soInd] += 1
                if 'char_info' in obs:
                    char_info = obs['char_info'][0]
                else:
                    # older (~2019) DRMs have char_status directly within obs
                    char_info = obs
                # must count planets for dets or chars, not all stars have both
                planet_num[soInd] = len(char_info['char_status'])
                # a characterization - any planet > 0 is an overall success
                # defs:
                #   [-1] success:  (>=1) planet is +1
                #   [-2] partial:  not success, and (>=1) planet is -1
                #   [-3] failure:  otherwise (i.e., char_status is all-0)
                success = np.any(char_info['char_status'] > 0)
                partial = ~success and np.any(char_info['char_status'] < 0)
                char_time = get_char_time(obs)
                if success:
                    obs_trace[soInd, nobs+1] = -1
                    char_ok.append((soInd, nobs))
                    char_success[soInd] = 1
                elif partial:
                    obs_trace[soInd, nobs+1] = -2
                    char_part.append((soInd, nobs))
                elif char_time == 0:
                    # a "miss" is an intTime=0 char from tieredScheduler
                    # it will look like a fail (status = 0) but have time = 0 too
                    char_miss.append((soInd, nobs))
                else:
                    obs_trace[soInd, nobs+1] = -3
                    char_nok.append((soInd, nobs))
            # Note: starshade-only missions have both char and det info for some observations
        # indexing
        self.Nsobs = Nsobs
        self.s2so = s2so
        # data reduction
        self.obs_trace = obs_trace
        self.det_ok    = det_ok
        self.det_nok   = det_nok
        self.det_nul   = det_nul
        self.det_fail  = det_fail
        self.det_iwa   = det_iwa
        self.det_owa   = det_owa
        self.char_ok   = char_ok
        self.char_nok  = char_nok
        self.char_part = char_part
        self.char_miss = char_miss
        self.arrival_times = arrival_times
        self.det_num_obs   = det_num_obs
        self.char_num_obs  = char_num_obs
        self.planet_num    = planet_num
        self.char_success  = char_success

########################################
###
###  Plot Container Class
###
########################################

class plotStarObsContainer(object):
    """A container class for methods relating to plotting star/observation traces
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

    def plot_star_obs_trace(self, sim):
        r'''Plot the star/obs trace.'''
        # one figure for whole plot
        # can try layout='constrained' for other placement rules
        fig = plt.figure(figsize=(18,20))
        gs = fig.add_gridspec(1, 3, wspace=0, width_ratios=(1,0.02,0.08))
        (ax, ax2, ax3) = gs.subplots(sharey=True)
        #ax = fig.add_subplot(111)
        self.plot_star_obs_core(fig, ax, ax2, ax3, sim)

        # overall title
        title = (
            'Observation Trace: Target Stars across Observations',
            'Detection Search Status: Blue Stripe; Characterization Status: Manila/Green Stripe',
            'Detections as Squares, Characterizations as Circles: See Legend')
        plt.suptitle('\n'.join(title), weight='bold', fontsize=20)
        # show the plot
        #plt.show(block=False)
        fname = 'star-obs-trace'
        plt.savefig(self.args.outpath % (fname, 'png'), **SAVEFIG_OPTS)
        # plt.savefig(self.args.outpath % (fname, 'pdf'))
        plt.close()


    def plot_star_obs_core(self, fig, ax, ax2, ax3, sim):
        # [0] -- colormap
        # allow for empty obs_trace
        obs_max = np.max(sim.obs_trace) if sim.obs_trace.size > 0 else 0.0
        vmax = max(8.0, obs_max) # let's say
        vmin = -3.0
        mapper = lambda x : (x - vmin)/(vmax - vmin)
        cmap = mpl.colors.LinearSegmentedColormap.from_list('Pooled', [
            (mapper(-3),   'lightcoral'),
            (mapper(-2),   'bisque'),
            (mapper(-1),   'lightgreen'),
            (mapper(0),    'white'),
            (mapper(vmax), 'royalblue')])
        mapper_plan = lambda x : (x - 0)/(2 - 0)
        cmap_plan = mpl.colors.LinearSegmentedColormap.from_list('Pooled', [
            (mapper_plan(0),   'white'),
            (mapper_plan(1),   'lightgreen'),
            (mapper_plan(2),   'forestgreen')])

        # [1] -- the base star-obs trace
        # print(f'uniques = {np.unique(sim.obs_trace)}')
        # star/obs trace serves as background for observations
        with warnings.catch_warnings():
            if sim.obs_trace.size == 0:
                warnings.simplefilter("ignore") # kill zero-size-image warning
            pcm = ax.imshow(sim.obs_trace, aspect='auto', origin='lower',
                            vmin=vmin, vmax=vmax, cmap=cmap)
        # (not working: steals space from ax, semi-fix with layout=constrained
        # in the figure() causes other issues with annotations)
        # fig.colorbar(pcm, ax=ax, shrink=0.6, location='bottom')

        # [2] -- overlay the individual observations
        # smaller markers, needed when Nsobs > 200 or so
        scalefac = 160.0 / max(160.0, sim.Nsobs)
        lw = 2.0*np.sqrt(scalefac); msD = 40*scalefac; msC = 80*scalefac
        obs_point_types = {
            'det_ok':    dict(label='D: Success',   c='g', marker='s', s=msD),
            'det_nul':   dict(label='D: No Planet', c='r', marker='s', s=msD),
            'det_iwa':   dict(label='D: Fail/IWA',  c='w', ec='r', lw=lw, marker='s', s=msD),
            'det_owa':   dict(label='D: Fail/OWA',  c='k', ec='r', lw=lw, marker='s', s=msD),
            'det_fail':  dict(label='D: Fail/SNR',  c='lime', ec='r', lw=lw, marker='s', s=msD),
            'char_ok':   dict(label='C: OK (Full)', c='forestgreen', marker='o', ec='k', s=msC),
            'char_part': dict(label='C: Partial',   c='orange', marker='o', ec='k', s=msC),
            'char_miss': dict(label='C: Miss (t=0)',c='lightgray', marker='o', ec='k', s=msC),
            'char_nok':  dict(label='C: Fail',      c='r', marker='o', ec='k', s=msC),
            }
        handles = []
        can_omit_point_types = ['char_miss']
        for pt_name, pt_props in obs_point_types.items():
            obs_list = getattr(sim, pt_name)
            if (not obs_list) and pt_name in can_omit_point_types:
                continue
            soInd_1 = [soInd for soInd, nobs in obs_list]
            nobs_1  = [nobs+0.5  for soInd, nobs in obs_list]
            h1 = ax.scatter(nobs_1, soInd_1, **pt_props)
            handles.append(h1)
        ax.legend(handles=handles)

        # [2b] -- side plot
        ax2.imshow(sim.planet_num.reshape((sim.Nsobs, 1)),  cmap=cmap_plan, 
                       origin='lower', aspect='auto', vmin=0, vmax=2)
        ax2.set_xticks([])
        ax2.set_title('#Planet', rotation=90, weight='bold', fontsize=14)
        planet_fail = np.logical_and((sim.char_success == 0), (sim.planet_num > 0))
        planet_fail_inx = np.nonzero(planet_fail)[0]
        # smaller markers, needed when Nsobs > 200 or so
        s_mark = int(80 * 70.0/max(70.0, sim.Nsobs))
        ax2.scatter(0*planet_fail_inx, planet_fail_inx, color='r', marker='x', s=s_mark)

        h_det  = ax3.barh(np.arange(sim.Nsobs), sim.det_num_obs, label='Det. Visit')
        h_char = ax3.barh(np.arange(sim.Nsobs), sim.char_num_obs, left=sim.det_num_obs, color='orange', label='Char. Visit')
        ax3.set_title('#Visit', rotation=90, weight='bold', fontsize=14)
        ax3.legend(handles=[h_det, h_char],  bbox_to_anchor=(1.04, 0.5), loc="center left")
        

        ## [3] -- decorations
        self.style_plot()

        if not sim.spc:
            y_lbl = [f'{sInd}' for sInd in sim.s2so.keys()]
        else:
            names = sim.spc['Name']
            y_lbl = [f'{names[sInd]}' for sInd in sim.s2so.keys()]

        # tick labels/locations
        Nstar = len(sim.s2so)
        # ensure y-axis tick labels are legible
        # (don't fail for Nstar == 0, which can happen)
        fs = min(14, max(4, int(1040/max(Nstar, 1))))
        pad = int(fs*6)
        ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))
        ax.set_yticks(range(sim.obs_trace.shape[0]), y_lbl, fontsize=fs, ha='left')
        ax.tick_params(pad=pad, axis='y') # left alignment => pad is needed (in points)
        # xticks
        x_top = sim.obs_trace.shape[1]
        ax.set_xlim(0, max(1, x_top-1))
        if sim.spc:
            xticks_given = ax.get_xticks()
            x_lbl = []
            x_tik = []
            for t in xticks_given:
                if t >= x_top:
                    continue
                try:
                    lbl = f'{t}\n{int(sim.arrival_times[int(t+1)])}d\n{sim.arrival_times[int(t+1)]/365.25:.2f}y'
                except IndexError:
                    lbl = f'{t}' # sometimes xticks_given[end] >> x_top
                x_tik.append(t)
                x_lbl.append(lbl)
            ax.set_xticks(x_tik, x_lbl, fontsize='large')

        # labels
        ax.set_ylabel('Target Name', weight='bold', fontsize=16)
        ax.set_xlabel('Observation Number or Mission Time', weight='bold', fontsize=16)
        if sim.show_seed:
            title = f'Mission Star-Observation Trace for {sim.name}'
        else:
            title = f'Mission Star-Observation Trace'
        ax.set_title(title, weight='bold', fontsize=18)


def main(args):
    # load the simulation run
    sim = SimulationRun(args.drm)
    # reduce DRM to get observation times
    sim.get_obs_trace()
    # open a plot container object
    plotter = plotStarObsContainer(args)
    # make plots
    plotter.plot_star_obs_trace(sim)
    print('%s: Done.' % args.progname)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Plot trace of target stars and observations.", epilog='')
    parser.add_argument('drm', metavar='DRM', default='', help='drm file')
    parser.add_argument('-o', default='./plot-%s.%s', type=str, dest='outpath', help='Output file pattern.')
    #parser.add_argument('-v', default=False, action='store_true', 
    #                        dest='verbose', help='Verbosity.')
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
    while not os.path.exists(directory) and tries < 2:
        try:
            os.makedirs(directory)
        except OSError:
            tries += 1
            time.sleep(random.uniform(0.1, 0.2))
            pass # hope it was just concurrent creation 

    # do it
    main(args)
    sys.exit(0)

