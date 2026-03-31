r'''
utils.py -- Common utilities for DRM reduction

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
    '''Load reduction config file from dirname as tool for RpLBins.

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
    if not os.access(fn, os.R_OK):
        return None
    # If it exists, it's an error for it to not load as a mapping
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
    return d
                

