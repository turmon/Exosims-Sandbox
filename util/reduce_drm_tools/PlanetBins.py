#!/usr/bin/env python
r'''
PlanetBins.py -- Class for binning planets by physical characteristics.

Centralized so binning can stay consistent across modules here.
    
Can run as an executable, to verify customization, like so:

```
$ PlanetBins.py [DIR]
```

where `DIR` is an optional directory to look for a `config-reduce.json`
customization file. It will print a summary of the bin parameters
after applying customization. A useful idiom is:

```
$ diff <(util/reduce_drm_tools/PlanetBins.py sims/SCENARIO) <(util/reduce_drm_tools/PlanetBins.py)
```

which will output any differences from the customization.
'''

# turmon mar 2026

import os
import sys
import json
import argparse
import textwrap
from pathlib import Path
from collections import defaultdict
import numpy as np
import astropy.units as u

# idiom for imports from "."
try:
    from . import utils
except ImportError:
    import utils
strip_units = utils.strip_units


class RpLBins:
    r'''Class to hold Rp - Luminosity bin properties.

    Both a container object for bin properties, and a method-holding class.'''
    # guard variable: set to true when customize_parameters() was called
    # ensures that class instances have received customization (which
    # can be {}), even if multiprocessing is in effect
    _customized = False

    # known group names and their allowed attribute names
    _custom_groups = {
        'earthlike': [
            'Earth_Rp_scaled',
            'Earth_SMA_lo',
            'Earth_SMA_hi',
            'Earth_Rp_hi',
            'Earth_Rp_lo'],
        'RpL_bins': [
            'Rp_bins',
            'L_bins']
        }

    # extracted from the original definition
    # TBD: scaling flag present, but interpretation unclear, see elsewhere in this file
    Earth_Rp_scaled = True # scale Earth_Rp_lo to obtain sloped lower boundary
    Earth_SMA_lo = 0.95
    Earth_SMA_hi = 1.67
    Earth_Rp_hi = 1.40
    Earth_Rp_lo = 0.80

    # valid keys for Rp/L binning
    # Bin the detected planets into types
    # 1: planet-radius bin-edges  [units = Earth radii]
    # Original (v1), 3 x 3 bins:
    #   Rp_bins = np.array([0.5, 1.4, 4.0, 14.3])
    # Later (v2, May 2018):
    # 5 x 3 bins, see Kopparapu et al, arxiv:1802.09602v1,
    # Table 1 and in particular Table 3 column 1, column 2 and Fig. 2:
    Rp_bins = np.array([0.5, 1.0, 1.75, 3.5, 6.0, 14.3])
    # 2: stellar luminosity bins, in hot -> cold order
    #    Luminosity and SMA are inter-convertible: SMA = 1/sqrt(L).
    #    NB: *decreasing ordering* in L <=> increasing in SMA
    #    So for example, for the outer bin boundary:
    #         L = .0030  [unit: Lsun]
    #      => SMA = 18.26  [unit: AU]
    #    Referencing the radius/SMA plot (rad-sma-rectangle-bin-plot):
    #       Rp_bins controls y-axis
    #       1/sqrt(L_bins) controls x-axis for each Rp_bins stripe
    # v1:
    # L_bins = np.array([
    #    [185, 1.5,  0.38, 0.0065],
    #    [185, 1.6,  0.42, 0.0065],
    #    [185, 1.55, 0.40, 0.0055]])
    # v2:
    L_bins = np.array([
        [182, 1.0,  0.28, 0.0035],
        [187, 1.12, 0.30, 0.0030],
        [188, 1.15, 0.32, 0.0030],
        [220, 1.65, 0.45, 0.0030],
        [220, 1.65, 0.40, 0.0025],
        ])
    ## TODO: set up on __init__ based on L_bins -- for now,
    ## we are restricted to 3x5
    # set up the bin-number map from (Rp_bin, L_bin) -> RpL_bin
    # this is a map from (int,int) -> int:
    #   yields 0 for input pairs outside the allowed range
    #        [namely, 1..len(Rp_bins) and 1..len(L_bins[i])]
    #   yields 1...9 otherwise, with 1, 2, 3 for the small planets (Rp bin number = 1).
    Rp_L_to_RpL_bin = defaultdict(int)
    # there are many ways to set this up: here is one.
    # (the below lines are enough for the old 3x3 setup,
    # and they work for the new 5x3 setup too)
    Rp_L_to_RpL_bin[(1,1)] = 1 # smallest radius, highest luminosity => L bin number 1
    Rp_L_to_RpL_bin[(1,2)] = 2
    Rp_L_to_RpL_bin[(1,3)] = 3
    Rp_L_to_RpL_bin[(2,1)] = 4
    Rp_L_to_RpL_bin[(2,2)] = 5
    Rp_L_to_RpL_bin[(2,3)] = 6
    Rp_L_to_RpL_bin[(3,1)] = 7
    Rp_L_to_RpL_bin[(3,2)] = 8
    Rp_L_to_RpL_bin[(3,3)] = 9
    # New setup has 5*3 bins due to two new radius bins: added below
    Rp_L_to_RpL_bin[(4,1)] = 10
    Rp_L_to_RpL_bin[(4,2)] = 11
    Rp_L_to_RpL_bin[(4,3)] = 12
    Rp_L_to_RpL_bin[(5,1)] = 13
    Rp_L_to_RpL_bin[(5,2)] = 14
    Rp_L_to_RpL_bin[(5,3)] = 15

    def __init__(self):
        # convenenience variables useful for plots
        if not self._customized:
            print(f'RpLBins: Error: Class was not customized.', file=sys.stderr)
            raise RuntimeError("RpLBins class used before customization")

        # parametric w/r/t Earth_Rp_lo_1AU
        # derived attributes: earthlike
        if self.Earth_Rp_scaled:
            # default is to scale lower Rp boundary to get "Nevada" shape
            # these quantities are here for the plotter to later access, so it does not
            # have to embed too much domain knowledge
            self.Earth_Rp_lo1 = self.Earth_Rp_lo/np.sqrt(self.Earth_SMA_lo) # left boundary
            self.Earth_Rp_lo2 = self.Earth_Rp_lo/np.sqrt(self.Earth_SMA_hi) # right boundary
        else:
            # TODO: We want to extend is_earthlike() to other planet classes, but to do so,
            # we need to determine the behavior w/r/t SMA we want in is_earthlike for these
            # other classes. Once we decide that, we can:
            #   -- modify is_earthlike to obtain the desired behavior
            #   -- generalize rad-sma-rectangle-bin-plot.py so that the
            #      correct boundary is drawn
            #   -- insert any other needed data here (that the plotter, etc., may want)
            #raise RuntimeError("Unimplemented scaling flag.")
            self.Earth_Rp_lo1 = self.Earth_Rp_lo # left boundary
            self.Earth_Rp_lo2 = self.Earth_Rp_lo # right boundary

        # derived attributes: RpL histogram
        # the below ":" selectors are correct for increasing ordering
        self.L_lo = self.L_bins[:,:-1]
        self.L_hi = self.L_bins[:,1:]
        # 1b: bin lo/hi edges, same size as the resulting histograms
        #     simple copies of entries of Rp_bins
        self.Rp_lo = np.outer(self.Rp_bins[:-1], np.ones((3,1))).ravel()
        self.Rp_hi = np.outer(self.Rp_bins[1:],  np.ones((3,1))).ravel()
        # total number of bins (e.g., 9 = 3*4-3, for a 3x3 histogram)
        RpL_bin_count = self.L_bins.size - (self.Rp_bins.size - 1)
        # radius/luminosity bin boundaries
        #   if there are 9 bins, there are 10 bin-edges, at 0.5, 1.5, ..., 9.5.
        #   this histogram drops the "0", or out-of-range, RpL region
        self.RpL_bin_edge_list = np.arange(0, RpL_bin_count+1) + 0.5

    
    @classmethod
    def customize_parameters(cls, mapping):
        '''Set up the class parameters using a custom mapping passed in.

        Return the list of unmatched attribute names (but only those within
        attribute groups we care about, like "earthlike".
        We assume the mapping is a dict.
        '''
        cls._customized = True
        fails = []
        # Look for attribute groups we care about
        for group_name, ok_attrs in cls._custom_groups.items():
            if group_name in mapping:
                # plug mapping in where possible, as attributes
                for attr_name, value in mapping[group_name].items():
                    if attr_name.startswith('_'):
                        continue # comment
                    elif attr_name not in ok_attrs:
                        fails.append(attr_name)
                    else:
                        #print(f"binner: Set {attr_name} to {value}")
                        # assumption: lists convert to np.array
                        if isinstance(value, list) and not isinstance(value, str):
                            setattr(cls, attr_name, np.array(value))
                        else:
                            setattr(cls, attr_name, value)
        # insist on 5 x 3 bins for now
        if 'RpL_bins' in mapping:
            assert cls.Rp_bins.shape == (5+1, ), f"Rp wrong size: {cls.Rp_bins.shape = }"
            assert cls.L_bins.shape == (5, 3+1), f"L wrong size: {cls.L_bins.shape = }"
        return fails


    def show(self):
        r'''Print the set-able attributes that we know about.'''
        for group_name, attrs in self._custom_groups.items():
            print(f'Group: {group_name}')
            for attr in attrs:
                print(f'  attribute: {attr}')
                print(textwrap.indent(
                    str(getattr(self, attr, 'uninitialized!')),
                    " "*4))

            
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

        # Note: this is an assumption, typically true for JPL scripts
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
        if self.Earth_Rp_scaled:
            # Default case: we scale the lower Rp boundary (only)
            #  -- this is the definition used by the SAG13 notion of "earthlike"
            #  -- this is the same logic in EXOSIMS SurveySimulation Prototype
            Rp_plan_lo_s = self.Earth_Rp_lo/np.sqrt(a_plan)
        else:
            # Special case: this special scaling is turned off
            #  -- this makes is_earthlike() have axis-parallel Rp/SMA boundaries
            #  -- thus allowing is_earthlike() to be turned to other purposes
            Rp_plan_lo_s = self.Earth_Rp_lo
        # We use the numpy versions so that plan_ind can be a numpy vector.
        return np.logical_and(
           np.logical_and(Rp_plan >= Rp_plan_lo_s, Rp_plan <= self.Earth_Rp_hi),
           np.logical_and(a_plan  >= self.Earth_SMA_lo, a_plan  <= self.Earth_SMA_hi))



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Summarize DRMs.",
                                     epilog='')
    parser.add_argument('dir', metavar='DIR', type=str, nargs='?', default='',
                            help='directory to initialize from')
    args = parser.parse_args()
    args.progname = os.path.basename(sys.argv[0])

    # if directory provided, load customization
    if args.dir:
        p = Path(args.dir)
        if not p.is_dir():
            print(f'{args.progname}: Directory {args.dir} is not a readable directory', file=sys.stderr)
            sys.exit(1)
        params = utils.load_reduce_config(p)
        if params is None:
            print(f'{args.progname}: No customization found! Empty customization used.')
            params = {}
        else:
            print(f'{args.progname}: Loaded customization from file.')
    else:
        params = {}
        print(f'{args.progname}: No customization given, using defaults.')

    RpLBins.customize_parameters(params)
    binner = RpLBins()
    print(f'{args.progname}: Binner summary:')
    binner.show()
    print(f'{args.progname}: Done.')
    sys.exit(0)


