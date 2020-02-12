#!/usr/bin/env python
r'''
locus-wa-dmag.py: provide locus of WA, dMag pairs of simulated universe

usage:
  locus-wa-dmag.py [ -q ] [ -p ] SCRIPT CSVFILE

where:
  SCRIPT is a json script to initialize Exosims
  CSVFILE is a CSV file (or file pattern) with input fields for the table
and:
  [warning: below is out of date]
  -p means plain text output versus the default markdown
  -q means to continue quietly if CSVFILE is not present

Output goes to stdout, which was a mistake.

Typical usage:
  PYTHONPATH=EXOSIMS:Local util/locus-wa-dmag.py Scripts/HabEx_4m_TSDDold_DD_TF17_promo6_20190123.exp/s_c1=0.17_c2=0.35_c6=0.04.json sims/HabEx_4m_TSDDold_DD_TF17_promo6_20190123.exp/s_c1=0.17_c2=0.35_c6=0.04/reduce-earth-char-list.csv

turmon mar 2019
'''

from __future__ import division
import argparse
import sys
import os
import csv
from collections import defaultdict, Counter
import numpy as np
import astropy.units as u
import EXOSIMS.MissionSim

# global verbosity mode, also usable for debugging print's
VERBOSITY = 0
# redirection
DEV_NULL = open(os.devnull, 'w')

########################################
###
###  Utility Functions
###
########################################

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

def load_csv(infile):
    # load the CSV, convert some types
    info = []
    try:
        with open(infile, 'rb') as csvfile:
            reader = csv.DictReader(csvfile, delimiter=',')
            for row in reader:
                info.append(row)
    except IOError:
        sys.stderr.write("Error: could not open `%s'.  Quitting." % infile)
        raise
    # (no conversion possible)
    if len(info) == 0: return info
    # find conversion (just string-to-number where possible)
    # (we just inspect the first row)
    c_func = {}
    for k in info[0].keys():
        try:
            float(info[0][k])
            c_func[k] = float
        except ValueError:
            c_func[k] = str
    # do the conversion
    new_info = []
    for index,row in enumerate(info):
        new_row = {}
        for k in row:
            new_row[k] = c_func[k](row[k])
        new_row['_index'] = index
        new_info.append(new_row)
    return new_info

########################################
###
###  Specialized Functions
###
########################################

def load_sim(script, seed):
    # e.g.:
    #  script = 'Scripts/HabEx_4m_TSDDold_DD_TF17_promo6_20190123.exp/s_c1=0.17_c2=0.35_c6=0.04.json'
    #  seed = 10185257
    xspecs = {'seed': seed, 'ensemble_mode': 'standalone'}
    with RedirectStreams(stdout=DEV_NULL):
        sim = EXOSIMS.MissionSim.MissionSim(script, **xspecs)
    return sim

def locus_sim_star(sim, sInd, pInds):
    # Kepler's law for period -- skip this for now
    # Ms = sim.TargetList.MsTrue[sInd]
    # mu = const.G*Ms
    # sma = np.amax(sim.SimulatedUniverse.a[sim_pInds])
    # T = 2.*np.pi*np.sqrt(sma**3/mu)

    # arbitrary parameters
    dt = 2.0*u.d
    Nt = 400
    # convenience
    SU = sim.SimulatedUniverse
    # planet indexes in the sim
    sim_pInds = np.where(SU.plan2star == sInd)[0]
    # within these, indexes of earth candidates
    earth_pInds = sim_pInds[pInds]
    # loop over times
    Np = len(pInds)
    WA   = np.zeros((Np, Nt))
    dMag = np.zeros((Np, Nt))
    for index in range(Nt):
        SU.propag_system(sInd, dt)
        # WA,dMag for the previously-chosen planets of the given star
        WA[  :,index] = SU.WA[  earth_pInds].to('mas')
        dMag[:,index] = SU.dMag[earth_pInds]
    return WA, dMag
        
def locus_sim(sim, star_list):
    r'''Locus of WA,dMag for all stars/planets and a single sim.'''
    # lists of matching length, containing at each entry:
    #    star, planets, WAs(t1..tN), dMags(t1..tN)
    sInds = []
    pInds = []
    WAs   = []
    dMags = []
    for sInd_1, pInd_1 in star_list.iteritems():
        # note, pInd_1 is a Counter
        WA_1, dMag_1 = locus_sim_star(sim, sInd_1, sorted(pInd_1.keys()))
        WAs.append(WA_1)
        dMags.append(dMag_1)
        pInds.append(pInd_1.keys())
        sInds.append([sInd_1]*len(pInd_1))
    # collapse each list-of-numbers to an overall table-of-numbers
    WA = np.concatenate(WAs)
    dMag = np.concatenate(dMags)
    sInd = np.concatenate(sInds)
    pInd = np.concatenate(pInds)
    return sInd, pInd, WA, dMag
        
# sim.TimeKeeping.currentTimeAbs
# sim.TimeKeeping.currentTimeNorm

def locus_all(args, planet_list):
    r'''Iterate over seeds: load sim, find planet positions, dump.'''
    is_first = True
    for seed, star_list in planet_list.iteritems():
        sys.stderr.write('Sim (%d): Load, ' % (seed, )); sys.stderr.flush()
        sim = load_sim(args.script, seed)
        sys.stderr.write('propagate, '); sys.stderr.flush()
        sInd, pInd, WA, dMag = locus_sim(sim, star_list)
        sys.stderr.write('dump.\n'); sys.stderr.flush()
        dump_locus(is_first, args.fp_out, seed, sInd, pInd, WA, dMag)
        is_first = False

def find_planets(char_table):
    r'''Scan the table of observations and extract map from seed->star->planets.'''
    # char_table is a list of dicts containing char information.
    # plist[seed][sInd][pInd] = number of times (seed, sInd, pInd) was characterized
    #   it is a defaultdict of defaultdict of Counter
    plist = defaultdict(lambda: defaultdict(Counter))
    for row in char_table:
        # filter out some rows
        if row['is_deep'] != 1: continue
        if row['is_success'] != 0: continue
        plist[
            int(row['ensemble'])][
                int(row['sind'])][int(row['pind'])] += 1
    return plist

def dump_locus(is_first, fp_out, seed, sInd, pInd, WA, dMag):
    r'''Dump the WA/dMag locus for one full sim.'''
    # header
    if is_first:
        Ntime = WA.shape[1]
        fp_out.write('ensemble,sInd,pInd,')
        fp_out.write(','.join('WA_%d' % (i+1) for i in range(Ntime)))
        fp_out.write(',')
        fp_out.write(','.join('dMag_%d' % (i+1) for i in range(Ntime)))
        fp_out.write('\n')
    for i in range(len(sInd)):
        fp_out.write('%d,' % seed)
        fp_out.write('%d,' % sInd[i])
        fp_out.write('%d,' % pInd[i])
        fp_out.write(','.join('%g' % x for x in WA[i,:]))
        fp_out.write(',')
        fp_out.write(','.join('%g' % x for x in dMag[i,:]))
        fp_out.write('\n')

########################################
###
###  Main routine
###
########################################

def main(args):
    # load CSV with all characterizations
    char_table = load_csv(args.info)
    # find ensembles, stars, planets observed
    planet_list = find_planets(char_table)
    #import pdb; pdb.set_trace()
    # find and dump the locus for each observation
    locus_all(args, planet_list)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Locus of planet positions.",
                                     epilog='')
    parser.add_argument('script', metavar='SCRIPT', help='script file name')
    parser.add_argument('info', metavar='INFO.csv', help='seed/star info file name')
    #parser.add_argument('-o', '--outfile', metavar='FILE', help='output file name (or template)',
    #                        default='')
    parser.add_argument('-v', help='verbosity', action='count', dest='verbose',
                            default=0)
    args = parser.parse_args()
    
    # set umask in hopes that files will be group-writable
    os.umask(0o002)

    VERBOSITY = args.verbose

    # program name, for convenience
    args.progname = os.path.basename(sys.argv[0])
    
    # output file
    args.fp_out = sys.stdout
    
    # is the named script there
    if not os.path.isfile(args.script):
        sys.stderr.write("%s: Script (%s) not readable.\n" % (
            args.progname, args.script))
        sys.exit(1)

    # is the named info there
    if not os.path.isfile(args.info):
        sys.stderr.write("%s: Info (%s) not readable.\n" % (
            args.progname, args.info))
        sys.exit(1)

    main(args)
    sys.exit(0)
