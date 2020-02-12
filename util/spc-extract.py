#!/usr/bin/env python
#
# spc-extract: Extract info from SPC files -- somewhat ad hoc
#
# For usage, use the -h option.  Some options may be described there but not here.
#
# Typical usage:
#   spc-extract.py sims/script/spc/*.spc
# where:
# optionally:
#  -o FILE -- output file


# history:
#  turmon jun 2019 - created


import sys
import glob
import argparse
import os
import csv
import cPickle as pickle
import numpy as np
import astropy.units as u
import astropy.constants as const


############################################################
#
# Utility Functions
#
############################################################

def strip_units(x):
    r'''Strip astropy units from x.'''
    # TODO: allow coercing units to a supplied value
    if hasattr(x, 'value'):
        return x.value
    else:
        return x

# Container class for loading canned star-planet configurations
class StarPlanetInfo(object):
    r"""Star-planet configuration, as loaded from an external pickle."""

    def is_earthlike_all(self):
        r'''Is the planet earthlike?'''
        # handy abbreviations
        spc = self.spc
        plan2star = spc['plan2star']
        # extract planet and star properties
        Rp_plan = strip_units(spc['Rp'])
        a_plan = strip_units(spc['a'])
        L_star = spc['L'][plan2star]
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
        Rp_plan_lo = 0.80/np.sqrt(a_plan)
        # We use the numpy versions so that plan_ind can be a numpy vector.
        return np.logical_and(
            np.logical_and(Rp_plan >= Rp_plan_lo, Rp_plan <= 1.4),
            np.logical_and(L_plan  >= 0.3586,     L_plan  <= 1.1080))


    def load_from_spc_file(self, spc):
        r'''Load the star/planet info from a "spc" file given as an argument.
        This spc file transfer is compatible the Exosims ipyparallel output.'''
        # these will fail noisily if there is no file present
        print 'Loading SPC from', spc
        # spc file contains a dict with many fields - save them all
        self.spc = pickle.load(open(spc))
        # seed extracted from filename
        self.seed = os.path.splitext(os.path.basename(spc))[0]

    def __init__(self, spc):
        # load DRM and Star-Planet info
        self.load_from_spc_file(spc)


def process(args, info):
    r'''Add some info to the existing SPC fields.'''
    spc = info.spc
    # determine the period of the planets
    mu = const.G*(spc['Mp'] + spc['MsTrue'][spc['plan2star']])
    T = 2.0 * np.pi * np.sqrt(spc['a']**3 / mu)
    info.spc['T'] = T.to('d')
    # determine effective luminosity
    L_star = spc['L'][spc['plan2star']]
    L_plan = L_star / ((spc['a'].value)**2) # adjust star luminosity by distance^2 in AU
    info.spc['Lp'] = L_plan
    # insert earthlike info - as an int, not a boolean
    info.spc['earth'] = info.is_earthlike_all().astype(int)
    

def dump(args, n, info):
    # dump some info
    if not args.outfile:
        return
    # first-file
    if n == 0:
        with open(args.outfile, 'w') as _:
            pass
    # fields to dump
    vector_fields = ['Rp', 'T', 'Mp', 'a', 'Lp', 'earth']
    extra_fields = ['seed']
    # note, appending
    with open(args.outfile, 'a') as csvfile:
        # OK to have extra fields in dict-for-row
        w = csv.DictWriter(csvfile, fieldnames=(extra_fields+vector_fields), extrasaction='ignore')
        if n == 0: w.writeheader()
        # make a dictionary mapping field -> value
        for i in range(len(info.spc['Rp'])):
            d = {key:strip_units(info.spc[key][i]) for key in vector_fields}
            d['seed'] = info.seed
            w.writerow(d)

############################################################
#
# Main routine
#
############################################################

def main(args):
    for n, fn in enumerate(args.spcs):
        info = StarPlanetInfo(fn)
        process(args, info)
        dump(args, n, info)

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Pull SPC info.")
    parser.add_argument('spcs', metavar='SPCs', nargs='*', help='SPC file(s)')
    parser.add_argument('-o', '--outfile', help='name of output file',
                      dest='outfile', metavar='FILE', default='spc.out')
    
    args = parser.parse_args()
    
    # set umask in hopes that files/dirs will be group-writable
    os.umask(0o002)

    main(args)


