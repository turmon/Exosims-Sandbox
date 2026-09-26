r'''
utils.py -- Common utilities for DRM reduction

'''

# turmon mar 2026

import os
import re
import sys
import json
import argparse
import textwrap
from pathlib import Path
from collections import defaultdict
import numpy as np
import astropy.units as u

# idiom for imports from ".", if needed
#try:
#    from . import AnotherReductionModule
#except ImportError:
#    import AnotherReductionModule

# data reduction, local config filename
REDUCTION_CONFIG = 'config-reduce.json'


def strip_units(x):
    r'''Strip astropy units from x.'''
    # TODO: allow coercing units to a supplied value
    if hasattr(x, 'value'):
        return x.value
    else:
        return x


# 
def load_reduce_config(dirname, log_origin=None):
    '''Load reduction config file from dirname or its parent as tool for RpLBins.

    dirname is a Path.
    
    Returns:
     - None if not present (not an error)
     - The dictionary, if it is present
     '''
    if log_origin is None:
        log_string = 'PlanetBins.py'
    else:
        log_string = log_origin

    fn = dirname / REDUCTION_CONFIG
    # OK for it not to exist
    if not (fn.is_file() and os.access(fn, os.R_OK)):
        dirname_alt = dirname.parent
        if dirname_alt.suffix in ('.fam', '.exp') and dirname_alt.is_dir():
            # look one level up
            fn = dirname_alt / REDUCTION_CONFIG
            if not (fn.is_file() and os.access(fn, os.R_OK)):
                return None
            else:
                pass # fn is readable -- continue
        else:
            return None
    # Below here: fn is a readable file, either in dirname or dirname.parent
            
    # If fn exists, it's an error for it to not load as a mapping
    try:
        with open(fn, 'r') as fp:
            d = json.load(fp)
    except FileNotFoundError:
        print(f'{log_origin}: Error: Reduction configuration ({fn}) exists but unreadable.')
        raise
    except json.JSONDecodeError:
        print(f'{log_origin}: Error: Could not read JSON in {fn}')
        raise
    if not isinstance(d, dict):
        print(f'{log_origin}: Error: Reduction configuration ({fn}) is not a mapping.')
        raise ValueError("Reduction configuration was not a mapping.")
    # record where it came from
    #    str() is the relative path starting from sims/
    d['_config_filename'] = str(fn)
    return d


# Spec-file (JSON script or outspec) lookup orders for infer_spec_for_drm.
# Names refer to these files, for drm = sims/ENS/drm/SEED.pkl:
#   outspec-seed:     sims/ENS/log/outspec/SEED.json  (per-seed outspec)
#   outspec-seed-old: sims/ENS/run/outspec_SEED.json  (per-seed outspec, older layout)
#   reduce-outspec:   sims/ENS/reduce-outspec.json    (outspec copied by reduce_drms.py)
#   reduce-script:    sims/ENS/reduce-script.json     (script copied by reduce_drms.py)
#   script:           Scripts/ENS.json                (the original script)
# The outspec has parameter values as actually used (including Exosims defaults),
# but it also pins machine-specific default paths (e.g., catalogpath, spkpath) from
# the run. So, to instantiate Exosims objects, the script is the safer choice.
# Note: util/plot-obs-timelines.sh has its own (bash) lookup, in the
# outspec-first order; keep the two consistent if either changes.
SPEC_ORDER_OUTSPEC_FIRST = ('outspec-seed', 'outspec-seed-old', 'reduce-outspec', 'reduce-script', 'script')
SPEC_ORDER_SCRIPT_FIRST  = ('reduce-script', 'script')


def spec_candidates(drm_path, order=SPEC_ORDER_SCRIPT_FIRST):
    r'''Return the list of candidate spec files (Paths) for a DRM, in the given order.

    See SPEC_ORDER_* above for the names used in order.'''
    drm = str(drm_path)
    # sims/ENS/drm/SEED.pkl -> sims/ENS, SEED
    simdir = Path(drm.split('/drm/')[0])
    seed = Path(drm).stem
    # sims/ENS/drm/SEED.pkl -> Scripts/ENS.json (as done in plot-obs-timelines.sh)
    script = Path(re.sub(r'.*sims/', 'Scripts/', drm.split('/drm/')[0] + '.json'))
    paths = {
        'outspec-seed':     simdir / 'log' / 'outspec' / f'{seed}.json',
        'outspec-seed-old': simdir / 'run' / f'outspec_{seed}.json',
        'reduce-outspec':   simdir / 'reduce-outspec.json',
        'reduce-script':    simdir / 'reduce-script.json',
        'script':           script,
        }
    return [paths[name] for name in order]


def infer_spec_for_drm(drm_path, order=SPEC_ORDER_SCRIPT_FIRST):
    r'''Find the spec file (JSON script or outspec) for a DRM, using Sandbox conventions.

    Returns (spec, tried): spec is the first readable candidate Path (None if
    there is none), and tried is the list of all candidate Paths, for messages.'''
    tried = spec_candidates(drm_path, order)
    for fn in tried:
        if fn.is_file() and os.access(fn, os.R_OK):
            return fn, tried
    return None, tried


