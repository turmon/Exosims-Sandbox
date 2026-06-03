#!/usr/bin/env python
"""spc-extract: Extract info from SPC files into CSV

For usage, use the -h option.  Some options may be described there but not here.

Typical usage:
  `spc-extract.py sims/script/spc/*.spc`
or:
  `spc-extract.py -k .len sims/script/spc/*.spc`

SPC files, directories of `.spc` files, or scenario/experiment directories
may be given; directories are expanded recursively to all `.spc` files within.

The `-k` option can be repeated to name particular SPC keys to output.
If no `-k` is given, `.default-star` is used.  Special KEY values:

+ `.len`            => give field names and their vector lengths (diagnostic)
+ `.name`           => give field names only (diagnostic)
+ `.all`            => output fields of the most-common non-scalar length,
                       plus all scalar fields.  Use `--like KEY` to target
                       a specific length (e.g., `--like Mp` for planet-length
                       fields, `--like L` for star-length fields).
+ `.default-star`   => output a standard set of star keys (default when no -k given)
+ `.default-planet` => output a standard set of planet keys

Identifier columns, prepended before data columns in order scenario, basename, seed:

+ `-N` => scenario name (e.g., `sims/coroSched_20231122` becomes `coroSched_20231122`)
+ `-B` => scenario basename (last path component of the scenario name)
+ `-s` => seed (numeric stem of the `.spc` filename)

Other options:

+ `--json`     => emit JSON array of objects instead of CSV
+ `--like KEY` => with `-k .all`, select fields of the same length as KEY
+ `-o FILE`    => output file (default: stdout)
"""


import sys
import glob
import argparse
import os
import csv
import json
from collections import defaultdict
import pickle
import numpy as np
import astropy.units as u
import astropy.constants as const


############################################################
#
# Utility Functions
#
############################################################

class NumpyEncoder(json.JSONEncoder):
    r"""Custom JSON encoder for numpy types."""
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, u.quantity.Quantity):
            return obj.value
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


def strip_units(x):
    r'''Strip astropy units from x.'''
    # TODO: allow coercing units to a supplied value
    if hasattr(x, 'value'):
        return x.value
    else:
        return x

def get_scenario(spc_path):
    r'''Get scenario name from SPC file path, honoring sandbox conventions.

    E.g., sims/coroSched_20231122/spc/161215293.spc -> coroSched_20231122
          sims/exp.fam/scenarioX/spc/42.spc         -> exp.fam/scenarioX'''
    if not spc_path.startswith('sims/'):
        reasonable = os.path.dirname(spc_path)
        if reasonable.endswith('/spc'):
            return reasonable[:-4]
        return reasonable
    f_tail = spc_path[5:]
    d = os.path.dirname(f_tail)
    if d.endswith('/spc'):
        return d[:-4]
    return f_tail


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
        self.spc = pickle.load(open(spc, 'rb'), encoding='latin1')
        # seed extracted from filename
        self.seed = os.path.splitext(os.path.basename(spc))[0]
        # filename
        self.filename = spc
        # scenario name and basename from path
        self.scenario = get_scenario(spc)
        self.basename = os.path.basename(self.scenario)

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
    

def expand_spc(spcs, progname):
    r'''Expand spc input args so that directories are descended into.'''
    def expand_dir(d):
        dx = []
        for root, dirs, files in os.walk(d):
            if root.endswith('/spc'):
                dx.extend(glob.glob(f'{root}/*.spc'))
                dirs[:] = []
                continue
            if 'spc' in dirs:
                dx.extend(glob.glob(f'{root}/spc/*.spc'))
                dirs[:] = []
                continue
            downs = [d for d in dirs if (
                d.endswith('.exp') or d.endswith('.fam') or
                os.path.isdir(f'{root}/{d}/spc'))]
            dirs[:] = downs
        return dx

    d_all = []
    for x in spcs:
        if os.path.isfile(x):
            d_all.append(x)
        elif os.path.isdir(x):
            d_all.extend(expand_dir(x))
        else:
            print(f'{progname}: Fatal. Could not access {x}.', file=sys.stderr)
            sys.exit(1)
    return d_all


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
    

def select_fields(spc, keys, target_len=None, warn_mismatch=False, progname=''):
    r'''Classify keys into vector fields and scalar fields for output.

    Returns (fields, fields_scalar).
    target_len: if given, vector fields are those whose length equals target_len.
    warn_mismatch: if True, warn when non-scalar keys have mixed lengths.'''
    key_by_len = defaultdict(list)
    for k in keys:
        if k not in spc:
            sys.stderr.write(f"{progname}: Could not find key='{k}', skipping.\n")
            continue
        key_by_len[get_length(spc[k])].append(k)

    scalars = key_by_len.get(1, [])
    non_scalar = {l: ks for l, ks in key_by_len.items() if l != 1}

    if target_len is not None:
        if target_len == 1:
            return [], scalars
        return non_scalar.get(target_len, []), scalars

    if not non_scalar:
        return [], scalars

    # find the length with the most keys; break ties by preferring larger length
    best_len = max(non_scalar, key=lambda l: (len(non_scalar[l]), l))

    if warn_mismatch and len(non_scalar) > 1:
        dropped = [k for l, ks in non_scalar.items() if l != best_len for k in ks]
        sys.stderr.write(
            f"{progname}: Warning: mixed vector lengths in requested keys; "
            f"dropping {dropped} (not length {best_len}).\n")

    return non_scalar[best_len], scalars


def dump_names(args, n, info):
    r'''Just write the field names to the output.

    Don't use CSV-writer, because there are no headings, and the field names
    may in principle vary across lines.'''
    # make the ordering invariant, esp. over calls
    keys = sorted(info.spc.keys())
    for ef in args.extra_fields:
        args.out_fp.write(f'{ef}={getattr(info, ef)}\n')
    for f in keys:
        xtra = f'[{str(info.spc[f].shape)}]' if isinstance(info.spc[f],np.ndarray) else ''
        args.out_fp.write(f'{f}{xtra}\n')
    

def dump_lengths(args, n, info):
    r'''Just write the field names to the output.

    Don't use CSV-writer, because there are no headings, and the field names
    may in principle vary across lines.'''
    # make the ordering invariant, esp. over calls
    keys = sorted(info.spc.keys())
    # note, appending to output file
    # header
    if n == 0:
        args.out_fp.write(','.join(args.extra_fields + list(keys)) + '\n')
    row = [str(getattr(info, ef)) for ef in args.extra_fields]
    row += [str(get_length(info.spc[f])) for f in keys]
    args.out_fp.write(','.join(row) + '\n')
    

def dump(args, n, info):
    # allow parse-only usage
    if not args.out_fp: return
    # support all keys
    if '.all' in args.key:
        if args.like:
            if args.like not in info.spc:
                sys.stderr.write(f"{args.progname}: Warning: --like key '{args.like}' not in SPC, ignoring.\n")
                target_len = None
            else:
                target_len = get_length(info.spc[args.like])
            fields, fields_scalar = select_fields(info.spc, info.spc.keys(),
                                                  target_len=target_len,
                                                  progname=args.progname)
        else:
            fields, fields_scalar = select_fields(info.spc, info.spc.keys(),
                                                  progname=args.progname)
    else:
        fields, fields_scalar = select_fields(info.spc, args.key, warn_mismatch=True,
                                              progname=args.progname)
    # extra identifier field(s) to dump
    extra_fields = args.extra_fields
    all_fields = extra_fields + fields + fields_scalar
    # set up CSV writer (JSON path defers output to main())
    if not args.json:
        w = csv.DictWriter(args.out_fp, fieldnames=all_fields, extrasaction='ignore')
        if n == 0: w.writeheader()
    # make a dictionary mapping field -> value
    num_entries = len(info.spc[fields[0]]) if fields else 1
    for i in range(num_entries):
        d  = {key:strip_units(info.spc[key][i]) for key in fields}
        d.update({key:strip_units(info.spc[key]) for key in fields_scalar})
        d['seed']     = info.seed
        d['scenario'] = info.scenario
        d['basename'] = info.basename
        if args.json:
            args.rows.append({f: d[f] for f in all_fields})
        else:
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

    args.rows = []
    open_output(args)
    for n, fn in enumerate(args.spcs):
        info = StarPlanetInfo(fn)
        process(args, info)
        dumper(args, n, info)
    if args.json and args.out_fp:
        json.dump(args.rows, args.out_fp, indent=2, cls=NumpyEncoder)

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Extract star-planet configuration (SPC) info and output as CSV.",
        epilog="Special -k KEY values: .all (most-common-length fields), .default-star (star keys, default when no -k given), .default-planet (planet keys), .name (field names only), .len (field name/length table)")
    parser.add_argument('spcs', metavar='SPC', nargs='*', help='SPC file(s)')
    parser.add_argument('-o', '--outfile', help='name of output file, default stdout',
                      dest='outfile', metavar='FILE', default='-')
    parser.add_argument('-k', '--key', action='append', default=[], help='repeat to get multiple keys', type=str)
    parser.add_argument('-s', '--seed',     action='store_true', default=False,
                        help='include seed in output')
    parser.add_argument('-N', '--name',     action='store_true', default=False,
                        help='include scenario name in output')
    parser.add_argument('-B', '--basename', action='store_true', default=False,
                        help='include scenario basename in output')
    parser.add_argument('--json', help='JSON output format', action='store_true',
                        dest='json', default=False)
    parser.add_argument('--like', metavar='KEY', default=None,
                        help='with -k .all, select attributes of the same length as KEY')
    
    args = parser.parse_args()
    args.progname = os.path.basename(sys.argv[0])

    # set umask in hopes that files/dirs will be group-writable
    os.umask(0o002)

    # special keys
    default_fields_star = ['Name', 'Spec', 'L', 'MsTrue', 'int_comp', 'dist']
    default_fields_planet = ['Rp', 'T', 'Mp', 'a', 'I', 'Lp', 'earth', 's', 'plan2star']
    if '.default-star' in args.key or len(args.key) == 0:
        args.key = [k for k in args.key if k != '.default-star']
        args.key += default_fields_star
    elif '.default-planet' in args.key:
        args.key = [k for k in args.key if k != '.default-planet']
        args.key += default_fields_planet

    args.extra_fields = (
        (['scenario'] if args.name     else []) +
        (['basename'] if args.basename else []) +
        (['seed']     if args.seed     else [])
    )
    args.spcs = expand_spc(args.spcs, args.progname)
    main(args)


