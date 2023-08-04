#!/usr/bin/env python
# 
# json-attr-edit.py: edit attributes of a JSON script
#
# Usage:
#    json-attr-edit.py [-f] [-b BASE] [-a REPL] [-s REPL] [-c DIR] [-o OUTPUT] SCRIPT ...
#
# where
#   -b BASE gives a basket of pre-defined transforms
#   -a REPL replaces of every instance of the attribute named in REPL
#   -s REPL replaces a specific (drilled-down) attribute named in REPL
#   -c DIR gives a new cache directory
#   -o OUTPUT gives an output filename pattern
#             (the script name is plugged into the single %s in OUTPUT)
#   -f forces the write of generated files, possibly over-writing existing files
#
# See -h for more.

# turmon 2023-jul-30

import sys
import json
from json import JSONEncoder
import os
import re
import ast
import argparse
from pathlib import Path
from numbers import Number
from collections.abc import Mapping
import numpy as np

############################################################
#
# Pre-defined baskets of filename transformations
#
# This code defines some static transformations, mapping
# property names to transformation-functions, for the "-b" option.

# root directories in various formats
GATTACA_ROOTDIR = '/scratch_lg/exo-yield/EXOSIMS_external_files/'
MUSTANG_ROOTDIR = '/proj/exep/rhonda/Sandbox/Parameters/EXOSIMS_external_files/'
VARIABLE_ROOTDIR = '$EXOSIMS_PARAMS/'

# conversion functions

def to_basename(k, v):
    # "EZ_distribution": "/scratch_lg/exo-yield/EXOSIMS_external_files/exoZodi/nominal_maxL_distribution-Dec2019.fits",
    return Path(v).name

# unwanted duplication here ... still experimenting
def sub_g2v_vanilla(k, v):
    return v.replace(GATTACA_ROOTDIR, VARIABLE_ROOTDIR)

def sub_m2v_vanilla(k, v):
    return v.replace(MUSTANG_ROOTDIR, VARIABLE_ROOTDIR)

def sub_g2m_vanilla(k, v):
    return v.replace(GATTACA_ROOTDIR, MUSTANG_ROOTDIR)

def sub_m2g_vanilla(k, v):
    return v.replace(MUSTANG_ROOTDIR, GATTACA_ROOTDIR)

# all filename-containing properties
SCRIPT_PATH_PROPS = [
    "wdsfilepath",
    "binaryleakfilepath",
    "occHIPs_no",
    "include_known_RV",
    "core_thruput",
    "occ_trans",
    "core_mean_intensity",
    "EZ_distribution",
    ]

# from: gattaca -> variable
XFORM_g2v = {prop: sub_g2v_vanilla for prop in SCRIPT_PATH_PROPS}
XFORM_g2v["cachedir"] = lambda k,v: "$HOME/.EXOSIMS/cache/"

# from: mustang -> variable
XFORM_m2v = {prop: sub_g2v_vanilla for prop in SCRIPT_PATH_PROPS}
XFORM_m2v["cachedir"] = lambda k,v: "$HOME/.EXOSIMS/cache/"

# from: gattaca -> mustang
XFORM_g2m = {prop: sub_g2m_vanilla for prop in SCRIPT_PATH_PROPS}
XFORM_g2m["cachedir"] = lambda k,v: "$HOME/.EXOSIMS/cache/"

# from: mustang -> gattaca
XFORM_m2g = {prop: sub_m2g_vanilla for prop in SCRIPT_PATH_PROPS}
XFORM_m2g["cachedir"] = lambda k,v: "/cache/"

# all known baskets of transforms
REGISTRY = {
    "null": dict(),
    "g2v": XFORM_g2v,
    "g2m": XFORM_g2m,
    "m2g": XFORM_m2g,
    }


############################################################
#
# Main application code

def apply_xforms(level, item, bases, sform, xform, verbose):
    r'''Recursively expand a dictionary following mappings sform, xform.

    Scan each component of the nested dictionary-and-list structure
    "item". If any dictionary (or sub-dictionary) key matches an entry
    in the mapping "d", we invoke the function d[key] on it.
    '''
    if isinstance(item, str):
        return item
    elif isinstance(item, Number):
        return item
    elif isinstance(item, list):
        # recursively apply sform/xform to each entry of the list
        if verbose: print(f'{"  "*level}  recurse[list]: {bases = }; {len(item) = }')
        return [apply_xforms(level+1, i, (*bases, str(index)), sform, xform, verbose)
                    for index, i in enumerate(item)]
    elif isinstance(item, Mapping):
        # apply sform/xform to any matching dictionary items
        val = dict()
        for k, v in item.items():
            kx = '.'.join((*bases, k))
            if verbose: print(f'{"  "*level}  sform key: {kx}')
            if kx in sform:
                if verbose: print(f'{"  "*level}  apply(s): {kx}')
                val[k] = sform[kx] # sform[] is a literal not a lambda
            elif k in xform:
                if verbose: print(f'{"  "*level}  apply(x): {k}')
                val[k] = xform[k](k, v)
            else:
                if verbose: print(f'{"  "*level}  recurse[dict]: {k}')
                val[k] = apply_xforms(level+1, v, (*bases, k), sform, xform, verbose)
        return val
    assert False, f'Cannot be reached: {item}'


def decode_literal_attrs(attrs, verbose):
    r'''decode a series of command-line attributes given in strings k=v'''
    # matches (whitespace)(identifier)(whitespace)[=:](value)
    # where:
    #    identifier can be missionLife or systems.0.QE, disallow -1
    #    a separator of = or :
    #    value is what is left
    p = re.compile('\s*([\w.]+)\s*([:=])(.*)')
    sform = dict()
    for attr in attrs:
        if verbose: print(f'  Decoding |{attr}|...')
        # kind message if not key:val format
        m = p.match(attr)
        if not m:
            print(f'Error converting {attr}, need "key=value" format', file=sys.stderr)
            raise
        k, sep, v = m.groups()
        if sep == ':':
            # a plain string: take v as-is
            val = v
        elif sep == '=':
            # try to eval v as an expression
            try:
                val = ast.literal_eval(v)
            except:
                print(f'Error interpreting {attr} with eval', file=sys.stderr)
                raise
        if verbose: print(f'  Got {k} -> {val} (of type: {type(val)})')
        sform[k] = val
    return sform


def make_xforms(args):
    if args.verbose:
        print(f'{args.progname}: Decoding arguments.')
    # 1: make xform dictionary
    xform_base = REGISTRY[args.base]
    if args.cachedir:
        xform_base['cachedir'] = lambda k,v: args.cachedir
    xform_args = decode_literal_attrs(args.attr_x, args.verbose)
    # make two-input lambda's from the literals just found
    xform_adds = {k:(lambda x1, x2: v) for k, v in xform_args.items()}
    # merge these two (xform_adds supersedes)
    xform = {**xform_base, **xform_adds}
    # 2: make sform dictionary
    sform = decode_literal_attrs(args.attr_s, args.verbose)
    return sform, xform


def dump(args, script, out_spec):
    if args.output in ('', '-'):
        stream = sys.stdout
    else:
        basename = os.path.basename(script)
        fn = args.output % basename
        if fn == script:
            print(f'{args.progname}: Error: Output equals input ({script}). Skipping.', file=sys.stderr)
            return 0
        if os.path.isfile(fn) and not args.force:
            print(f'{args.progname}: Error: Output file exists ({fn}). Skipping; consider -f.', file=sys.stderr)
            return 0
        fn_dir = os.path.dirname(fn) # = '' if no dir/ present
        if fn_dir and not os.path.isdir(fn_dir):
            os.makedirs(fn_dir)
        stream = open(fn, 'w')
    json.dump(out_spec, stream, indent=2)
    if stream is not sys.stdout:
        stream.close()
    return 1


def main(args):
    all_ok = True
    try:
        sform, xform = make_xforms(args)
    except:
        print(f'Error converting input arguments, consider -v', file=sys.stderr)
        raise
    for script in args.script:
        if args.verbose:
            print(f'{args.progname}: Processing script = {script}')
        # get input script
        in_spec = json.load(open(script, 'r'))
        # transform the specs dictionary
        out_spec = apply_xforms(0, in_spec, (), sform, xform, args.verbose)
        # write it out
        ok = dump(args, script, out_spec)
        if not ok: all_ok = False
    return all_ok

############################################################

REPL_NOTE = '''
BASE is a basket of pre-selected transforms, default is none. Otherwise:
  * g2v is gattaca filename base to $EXOSIMS_PARAMS variable
  * g2m is gattaca filename base to mustang filename base
  * m2g is mustang filename base to gattaca filename base

REPL may be given in two ways:
  * a string of the form ATTR:VALUE, meaning that the named ATTR will
    be assigned the value VALUE as a string type (i.e., quoted in the JSON).
    e.g.:  modules.StarCatalog:EXOCAT1 sets "StarCatalog" to "EXOCAT1"
  * a string of the form ATTR=VALUE, meaning that the named ATTR will
    be assigned the value VALUE, where the string VALUE is eval'ed by
    Python, thus having a floating-point, integer, boolean, list, or other type.
    e.g.:  pupilDiam=6.0 sets "pupilDiam" to the number, 6.0
           lucky_planets=True sets "lucky_planets" to the boolean, True
           Rprange=[0.5,4.0] sets "Rprange" to that numeric vector
The above avoids most nested quoting by signaling strings with ":".

The above ATTR may be a bare name (eta, scaleOrbits, etc.), or a dotted name,
where dots stand for references. So, scienceInstruments.1.QE=0.72 means the second
scienceInstrument's QE is set to the floating-point value 0.72.
'''

############################################################

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Rename entries of JSON scripts.",
                                         epilog=REPL_NOTE,
                                         formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('script', metavar='SCRIPT', nargs='*', default=[], help='JSON script')
    parser.add_argument('-b', '--base', type=str, default='null', help=f'transform base, one of: {"|".join(REGISTRY.keys())}')
    parser.add_argument('-a', '--attr', type=str, action='append', default=[], metavar='REPL',
                           dest='attr_x', help='named attribute within script (e.g., missionTime)')
    parser.add_argument('-s', '--specific', type=str, action='append', default=[], metavar='REPL',
                           dest='attr_s', help='specific attribute within script (e.g., scienceInstruments.1.QE)')
    parser.add_argument('-c', '--cachedir', type=str, default='', metavar='DIR',
                            help='"cachedir" attribute of script (like -a cachedir)')
    parser.add_argument('-o', '--output', default='./xform-%s', type=str, help='Output file pattern (contains one %%s)')
    parser.add_argument('-f', '--force', default=False, action='store_true', help='Allow outputs to overwrite.')
    parser.add_argument('-v', default=False, action='store_true', dest='verbose', help='Verbosity.')
    args = parser.parse_args()
    args.progname = os.path.basename(sys.argv[0])
    
    if args.output not in ('', '-'):
        if args.output.count('%s') != 1:
            sys.stderr.write(f'{args.progname}: Fatal: Output {args.output} must have one %s pattern.')
            sys.exit(1)

    # set umask in hopes that files will be group-writable
    os.umask(0o002)

    # do it
    sys.exit(0 if main(args) == True else 1)
