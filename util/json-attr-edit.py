#!/usr/bin/env python
# 
# json-attr-edit.py: edit attributes of a JSON script
#
# Usage:
#    json-attr-edit.py [-f] [-v] [-b BASE] [-a REPL] [-s REPL] [-c DIR] [-o OUTPUT] SCRIPT ...
#
# where
#   -h gives more complete documentation
#   -b BASE gives a basket of pre-defined transforms (see -h for all)
#   -a REPL replaces of every instance of the attribute named in REPL
#   -s REPL replaces a specific (drilled-down) attribute named in REPL
#   -c DIR gives a new cache directory
#   -o OUTPUT gives an output filename pattern
#             (the script name is plugged into the single %s in OUTPUT)
#   -f forces the write of generated files, possibly over-writing existing files
#   -v increases verbosity
#
# See -h for more.

# turmon 2023-jul-30

import sys
import json
import os
import re
import ast
import argparse
import copy
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
#GATTACA_ROOTDIR = 'EXOSIMS_external_files/'
MUSTANG_ROOTDIR = '/proj/exep/rhonda/Sandbox/Parameters/EXOSIMS_external_files/'
VARIABLE_ROOTDIR = '$EXOSIMS_PARAMS/'

# old mustang conventions (say, before 2022?)
OLDMUST_ROOTDIR = '/proj/exep/rhonda/'
OLDMUST_ROOTOPT = '/proj/exep/rhonda/HabEx/'

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

# for this one (old mustang conventions -> variable), we have
# to special-case the starshade vs. coronagraph
# the names here are like:
#  "binaryleakfilepath": "/proj/exep/rhonda/binary_leakage/leakData.csv",
#  "wdsfilepath": "/proj/exep/rhonda/binary_leakage/ExoCat_WDS_Sep_dM.csv",
#  "EZ_distribution": "/proj/exep/rhonda/Sandbox/HabEx/EXOSIMS/EXOSIMS/ZodiacalLight/nominal_maxL_distribution-Dec2019.fits",
#  "occHIPs_no":"/proj/exep/rhonda/topStars/DD8_HabEx_4m_20190725.txt",  
#  "include_known_RV": "/proj/exep/rhonda/topStars/EPRV_stars/HIP_RV_musts.txt",
#  "occ_trans": "/proj/exep/rhonda/HabEx/Krist_occ_trans2_6m_asec500nm_OWA128.fits",
#  "core_thruput": "/proj/exep/rhonda/HabEx/Krist_core_thruput_6m_asec500nm_OWA128.fits", 
#  "occ_trans": "/proj/exep/rhonda/HabEx/Starshade/TV3_occ_trans_asec_650_60m_95200_IWA65.fits",
#  "core_thruput": "/proj/exep/rhonda/HabEx/Starshade/TV3_core_thruput_asec_650_60m_95200_IWA65.fits",
#  "core_mean_intensity": "/proj/exep/rhonda/HabEx/Starshade/TV3_core_mean_intensity....fits",

def sub_om2v_optics(k,v):
    if k == 'EZ_distribution':
        return VARIABLE_ROOTDIR + 'exoZodi/' + to_basename(k, v)
    elif '/HabEx/' in v:
        # an optical parameter:
        #   occ_trans, core_thruput, core_mean_intensity
        if 'Starshade' in v:
            return v.replace(OLDMUST_ROOTOPT + 'Starshade/', VARIABLE_ROOTDIR + 'starshade/')
        else:
            return v.replace(OLDMUST_ROOTOPT, VARIABLE_ROOTDIR + 'coronagraph/')
    else:
        return v.replace(OLDMUST_ROOTDIR, VARIABLE_ROOTDIR)

## enter a new basket:
##  - define a conversion subroutine for each keyword that will be altered (above)
##  - place them into a dictionary indexed by key name
##  - add a #note field to the dictionary
##  - enter in the REGISTRY variable below

# from: gattaca -> variable
XFORM_g2v = {prop: sub_g2v_vanilla for prop in SCRIPT_PATH_PROPS}
XFORM_g2v['#note'] = 'Gattaca 2023 filename convention -> EXOSIMS_PARAMS'
XFORM_g2v["cachedir"] = lambda k,v: "$HOME/.EXOSIMS/cache/"

# from: mustang -> variable
XFORM_m2v = {prop: sub_m2v_vanilla for prop in SCRIPT_PATH_PROPS}
XFORM_m2v['#note'] = 'Mustang 2023 filename convention -> EXOSIMS_PARAMS'
XFORM_m2v["cachedir"] = lambda k,v: "$HOME/.EXOSIMS/cache/"

# from: gattaca -> mustang
XFORM_g2m = {prop: sub_g2m_vanilla for prop in SCRIPT_PATH_PROPS}
XFORM_g2m['#note'] = 'Gattaca 2023 filename convention -> mustang'
XFORM_g2m["cachedir"] = lambda k,v: "$HOME/.EXOSIMS/cache/"

# from: mustang -> gattaca
XFORM_m2g = {prop: sub_m2g_vanilla for prop in SCRIPT_PATH_PROPS}
XFORM_m2g['#note'] = 'Mustang 2023 filename convention -> gattaca'
XFORM_m2g["cachedir"] = lambda k,v: "/cache/"

# from: old_mustang -> variable
XFORM_om2v = {prop: sub_om2v_optics for prop in SCRIPT_PATH_PROPS}
XFORM_om2v['#note'] = 'Mustang pre-2023 filename convention -> EXOSIMS_PARAMS'
XFORM_om2v["cachedir"] = lambda k,v: "$HOME/.EXOSIMS/cache/"

# all known baskets of transforms
REGISTRY = {
    "null": dict(),
    "g2v":  XFORM_g2v,
    "g2m":  XFORM_g2m,
    "m2g":  XFORM_m2g,
    "m2v":  XFORM_m2v,
    "om2v": XFORM_om2v,
    }


############################################################
#
# Main application code

def apply_xforms(level, item, bases, sform, eform, xform, verbose):
    r'''Recursively expand a dictionary following mappings sform, xform.

    Scan each component of the nested dictionary-and-list structure
    "item". If any dictionary (or sub-dictionary) key matches an entry
    in the mapping "d", we invoke the function d[key] on it.
    '''
    global SPECS
    if isinstance(item, str):
        return item
    elif isinstance(item, Number):
        return item
    elif isinstance(item, list):
        # recursively apply sform/xform to each entry of the list
        if verbose > 1: print(f'{"  "*level}  recurse[list]: {bases = }; {len(item) = }')
        return [
            apply_xforms(level+1, i, (*bases, str(index)), sform, eform, xform, verbose)
            for index, i in enumerate(item)]
    elif isinstance(item, Mapping):
        # apply sform/xform to any matching dictionary items
        val = dict()
        for k, v in item.items():
            # extended key name (scienceInstruments.0.QE)
            kx = '.'.join((*bases, k))
            if verbose > 1: print(f'{"  "*level}  sform key: {kx}')
            if kx in sform:
                if verbose > 1: print(f'{"  "*level}  apply(s): {kx}')
                val[k] = sform[kx] # sform[] is a literal not a lambda
                if verbose > 0: print(f'{kx}: {val[k]}')
            elif kx in eform:
                val[k] = eval(eform[kx], {"np": np, **SPECS})
                if verbose > 0: print(f'{kx}: {val[k]}')
            elif k in xform:
                if verbose > 1: print(f'{"  "*level}  apply(x): {k}')
                val[k] = xform[k](k, v)
                if verbose > 0: print(f'{k}: {val[k]}')
            else:
                if verbose > 1: print(f'{"  "*level}  recurse[dict]: {k}')
                val[k] = apply_xforms(level+1, v, (*bases, k), sform, eform, xform, verbose)
        return val
    assert False, f'Cannot be reached: {item}'


def decode_literal_attrs(attrs, verbose, force_string=False):
    r'''decode a series of command-line attributes given in strings k=v'''
    # matches (whitespace)(identifier)(whitespace)[=:](value)
    # where:
    #    identifier can be missionLife or systems.0.QE, disallow -1
    #    a separator of = or :
    #    value is what is left
    p = re.compile('\s*([\w.]+)\s*([:=])(.*)')
    sform = dict()
    for attr in attrs:
        if verbose > 1: print(f'  Decoding |{attr}|...')
        # kind message if not key:val format
        m = p.match(attr)
        if not m:
            print(f'Error converting {attr}, need "key=value" format', file=sys.stderr)
            raise
        k, sep, v = m.groups()
        if sep == ':' or force_string:
            # a plain string: take v as-is
            val = v
        elif sep == '=':
            # try to eval v as an expression
            try:
                val = ast.literal_eval(v)
            except:
                print(f'Error interpreting {attr} with eval', file=sys.stderr)
                raise
        if verbose > 1: print(f'  Got {k} -> {val} (of type: {type(val)})')
        sform[k] = val
    return sform


def make_xforms(args):
    if args.verbose > 1:
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
    # 2: make eform dictionary
    eform = decode_literal_attrs(args.attr_e, args.verbose, force_string=True)
    return sform, eform, xform


def dump(args, script, out_spec):
    if args.output in ('', '-'):
        stream = sys.stdout
    else:
        basename = os.path.basename(script)
        if '%s' in args.output:
            fn = args.output % basename
        else:
            fn = args.output.format(**{'__script__': basename, **out_spec})
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
    global SPECS
    all_ok = True
    try:
        sform, eform, xform = make_xforms(args)
    except:
        print(f'Error converting input arguments, consider -v', file=sys.stderr)
        raise
    for script in args.script:
        if args.verbose > 0:
            print(f'{args.progname}: Processing script = {script}')
        # get input script
        in_spec = json.load(open(script, 'r'))
        SPECS = copy.deepcopy(in_spec)
        # transform the specs dictionary
        out_spec = apply_xforms(0, in_spec, (), sform, eform, xform, args.verbose)
        # write it out
        ok = dump(args, script, out_spec)
        if not ok: all_ok = False
    return all_ok

############################################################

REPL_NOTE_template = '''
BASE is an optional basket of pre-selected transforms, one of:
  * |DOCS|

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
The first way avoids nested quoting by signaling strings with ":".

The above ATTR may be a bare name (eta, scaleOrbits, etc.), or a dotted name,
where dots stand for references. So, scienceInstruments.1.QE=0.72 means the second
scienceInstrument's QE is set to the floating-point value 0.72.

Supply "-v" to show substitutions that were made.
'''

# substitute in the (dynamic) replacement bases
REPL_NOTE = REPL_NOTE_template.replace(
    '|DOCS|', 
    '\n  * '.join(f'{name}: {base["#note"]}' for name,base in REGISTRY.items() if base))

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
    parser.add_argument('-e', '--eval', type=str, action='append', default=[], metavar='REPL',
                           dest='attr_e', help='specific attribute within script (e.g., scienceInstruments.1.QE)')
    parser.add_argument('-c', '--cachedir', type=str, default='', metavar='DIR',
                            help='"cachedir" attribute of script (like -a cachedir:DIR)')
    parser.add_argument('-o', '--output', default='./xform-%s', type=str, help='Output file pattern (contains one %%s)')
    parser.add_argument('-f', '--force', default=False, action='store_true', help='Allow outputs to overwrite.')
    parser.add_argument('-v', default=0, action='count', dest='verbose', help='Verbosity.')
    args = parser.parse_args()
    args.progname = os.path.basename(sys.argv[0])
    
    if args.output not in ('', '-'):
        if args.output.count('%s') != 1 and '{' not in args.output:
            sys.stderr.write(f'{args.progname}: Fatal: Output {args.output} must have one formatting pattern.')
            sys.exit(1)

    # set umask in hopes that files will be group-writable
    os.umask(0o002)

    # do it
    sys.exit(0 if main(args) == True else 1)
