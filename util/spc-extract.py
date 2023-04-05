#!/usr/bin/env python
#
# spc-extract: Extract info from SPC files into CSV
#
# For usage, use the -h option.  Some options may be described there but not here.
#
# Typical usage:
#   spc-extract.py sims/script/spc/*.spc
# or:
#   spc-extract.py -k .len sims/script/spc/*.spc
#
# where -k can be used more than once to name particular keys to output,
# or can take a handful of special values:
#   .len  => give key names and the corresponding vector lengths
#   .name => give just key names
#   .all  => output as many keys as possible in one CSV
#            (outputs the keys that correspond to the most-commonly-seen
#             vector length, e.g., stars)
#   .next => output the keys corresponding to the second-most-seen vector length
#            (e.g., planets).  Not implemented.
#   .default-star, .default-planet => output typically-present keys for stars
#            and planets, respectively
#
# optionally:
#  -o FILE -- output file


# history:
#  turmon jun 2019 - created
#  turmon sep 2020 - refined


from __future__ import print_function
import sys
import glob
import argparse
import os
import csv
from collections import defaultdict
import six.moves.cPickle as pickle
import numpy as np
import astropy.units as u
import astropy.constants as const
from six.moves import range

# unpickling python2/numpy pickles within python3 requires this
PICKLE_ARGS = {} if sys.version_info.major < 3 else {'encoding': 'latin1'}


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
        if False:
            print('Loading SPC from', spc)
        # spc file contains a dict with many fields - save them all
        self.spc = pickle.load(open(spc, 'rb'), **PICKLE_ARGS)
        # seed extracted from filename
        self.seed = os.path.splitext(os.path.basename(spc))[0]
        # filename
        self.filename = spc
        self.basename = os.path.basename(spc)

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
    # adjust star luminosity by distance^2 in AU
    L_plan = L_star / ((spc['a'].value)**2)
    info.spc['Lp'] = L_plan
    # insert earthlike info - as an int, not a boolean
    info.spc['earth'] = info.is_earthlike_all().astype(int)
    

def open_output(args):
    r'''Prepare the output file.'''
    if not args.outfile:
        # skip it
        args.out_fp = None
    elif args.outfile == '-':
        # use stdout
        args.out_fp = sys.stdout
    else:
        # truncate output file
        args.out_fp = open(args.outfile, 'w')


def get_length(qty):
    r'''Get length of a certain quantity, 1 if scalar.'''
    try:
        l = len(qty)
    except TypeError:
        l = 1
    # not OK for bytes, maybe that's good
    if isinstance(qty, str):
        l = 1
    return l
    

def majority_key(spc, keys):
    r'''Get the keys sharing the most common field length.

    E.g., if T, Mp, and Rp have a common length, and two other fields
    have length 1, return these keys.'''
    # 
    key_by_len = defaultdict(list)
    for k in keys:
        if k not in spc:
            sys.stderr.write("Could not find key=`%s', skipping.\n" % k)
            continue
        key_by_len[get_length(spc[k])].append(k)
    # find most common length
    l_max = 0
    for l, klist in key_by_len.items():
        if len(klist) > l_max:
            l_max = l
        # in general, longer vector wins ties
        #    and in particular, vector keys win ties with scalars
        if len(klist) == l_max and l > l_max:
            l_max = l
    # get the keys...determine the fields for scalars, too
    if l_max == 1:
        fields = []
        fields_scalar = key_by_len[l_max]
    else:
        fields = key_by_len[l_max]
        fields_scalar = key_by_len[1]
    # return both
    return fields, fields_scalar


def dump_names(args, n, info):
    r'''Just write the field names to the output.

    Don't use CSV-writer, because there are no headings, and the field names
    may in principle vary across lines.'''
    # make the ordering invariant, esp. over calls
    keys = sorted(info.spc.keys())
    args.out_fp.write('%s' % info.seed)
    for f in keys:
        args.out_fp.write(',%s' % f)
    args.out_fp.write('\n')
    

def dump_lengths(args, n, info):
    r'''Just write the field names to the output.

    Don't use CSV-writer, because there are no headings, and the field names
    may in principle vary across lines.'''
    # make the ordering invariant, esp. over calls
    keys = sorted(info.spc.keys())
    # note, appending to output file
    # header
    if n == 0:
        args.out_fp.write('seed')
        for f in keys:
            args.out_fp.write(',%s' % f)
        args.out_fp.write('\n')
    args.out_fp.write('%s' % info.seed)
    for f in keys:
        args.out_fp.write(',%s' % get_length(info.spc[f]))
    args.out_fp.write('\n')
    

def dump(args, n, info):
    # allow parse-only usage
    if not args.out_fp: return
    # support all keys
    if '.all' in args.key:
        fields, fields_scalar = majority_key(info.spc, info.spc.keys())
    else:
        # fields, fields_scalar = args.key, []
        fields, fields_scalar = majority_key(info.spc, args.key)
    # extra field(s) to dump
    extra_fields = ['seed']
    # (note, appending to output file)
    # OK to have extra fields in dict-for-row
    w = csv.DictWriter(args.out_fp, fieldnames=(extra_fields+fields+fields_scalar), extrasaction='ignore')
    if n == 0: w.writeheader()
    # make a dictionary mapping field -> value
    num_entries = len(info.spc[fields[0]]) if fields else 1
    for i in range(num_entries):
        d  = {key:strip_units(info.spc[key][i]) for key in fields}
        d.update({key:strip_units(info.spc[key]) for key in fields_scalar})
        d['seed'] = info.seed
        w.writerow(d)


############################################################
#
# Main routine
#
############################################################

def main(args):
    if '.name' in args.key:
        dumper = dump_names
    elif '.len' in args.key:
        dumper = dump_lengths
    else:
        dumper = dump

    open_output(args)
    for n, fn in enumerate(args.spcs):
        info = StarPlanetInfo(fn)
        process(args, info)
        dumper(args, n, info)

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Extract star-planet configuration (SPC) info and output as CSV.",
        epilog="Use KEY of .all for as many fields as possible, .name for fieldnames only, .len for field lengths")
    parser.add_argument('spcs', metavar='SPC', nargs='*', help='SPC file(s)')
    parser.add_argument('-o', '--outfile', help='name of output file, default stdout',
                      dest='outfile', metavar='FILE', default='-')
    parser.add_argument('-k', '--key', action='append', default=[], help='repeat to get multiple keys', type=str)
    
    args = parser.parse_args()
    
    # set umask in hopes that files/dirs will be group-writable
    os.umask(0o002)

    # special keys
    default_fields_star = ['Name', 'Spec', 'L', 'MsTrue', 'comp0', 'dist', 's']
    default_fields_planet = ['Rp', 'T', 'Mp', 'a', 'I', 'Lp', 'earth', 'plan2star']
    if '.default-star' in args.key or len(args.key) == 0:
        args.key = [k for k in args.key if k != '.default-star']
        args.key += default_fields_star
    elif '.default-planet' in args.key:
        args.key = [k for k in args.key if k != '.default-planet']
        args.key += default_fields_planet

    main(args)


