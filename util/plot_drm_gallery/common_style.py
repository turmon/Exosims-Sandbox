#!/usr/bin/env python
"""
Common styles, functions, etc.

"""

import os
import sys
import pandas as pd

# reduce_drm_tools lives one level up, in util/.  Make it importable whether we
# were started by the driver (which puts util/ on sys.path) or standalone from
# within this directory.
_UTIL_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _UTIL_DIR not in sys.path:
    sys.path.append(_UTIL_DIR)
from reduce_drm_tools.PlanetNames import PlanetNames


class PlotTracker:
    """Track graphics files written by a plot routine."""

    def __init__(self, ext_list=None):
        # list of filenames written
        self._files = []               
        # file extensions to make
        if ext_list is None:
            self._ext_list = ['png'] # default
        else:
            self._ext_list = ext_list[:]

    def set_ext_list(self, ext_list):
        """Set the file extensions for subsequent writes."""
        self._ext_list = list(ext_list)

    def write_plots(self, fig, dest_name, dest_tmpl,
                    ext_list=None, verbose=True, dpi=200, facecolor='none'):
        """Write figure to file(s) and record what was written."""
        if ext_list is None:
            ext_list = self._ext_list
        for ext in ext_list:
            fn_gfx = dest_tmpl % (dest_name, ext)
            if verbose:
                print(f'\tExport: {os.path.basename(fn_gfx)}')
            if facecolor is not None:
                fig.patch.set_facecolor(facecolor)
            fig.savefig(fn_gfx, dpi=dpi, bbox_inches='tight')
            self._files.append(os.path.basename(fn_gfx))

    def get_files(self):
        """Return list of filenames written."""
        return list(self._files)


# Helper function to create title from reduce_info dict
def plot_make_title(reduce_info):
    """Create plot title from metadata dict"""
    if not reduce_info:
        rv = ''
    elif 'experiment' in reduce_info:
        exp_name = str.strip(reduce_info['experiment'])
        if len(exp_name) < 50:
            chaser = f", Ensemble Size {reduce_info['ensemble_size']}"
        else:
            chaser = ''
        rv = exp_name + chaser
    else:
        rv = ''
    return rv


def load_csv_files(src_tmpl, csv_files):
    """Load CSV files and return as a list of DataFrames.

    For use by standalone plot scripts (not the driver). Exits on failure.
    """
    dataframes = []
    for csv_name in csv_files:
        csv_path = src_tmpl % (csv_name, 'csv')
        try:
            dataframes.append(pd.read_csv(csv_path))
        except Exception as e:
            print(f"Fatal: Could not load CSV file '{csv_path}': {e}", file=sys.stderr)
            sys.exit(1)
    return dataframes



def planet_names(reduce_info):
    """Return the PlanetNames (display names of the earthlike planet class).

    The class can be re-defined numerically in config-reduce.json, in which
    case calling it an "Earth" in a label is wrong.  Falls back to the
    historical Earth-based names when nothing was customized.
    """
    return PlanetNames.from_reduce_info(reduce_info)


def load_reduce_info(src_tmpl):
    """Load reduce-info.csv as a dict, plus the planet-class display names.

    The names come from the config-reduce.json reachable from the scenario
    directory implied by src_tmpl.  Used both by plot_drm_driver.py and by
    each module's standalone main(), so the two entry points agree.

    Note: the names ride along as plain strings, so they survive pickling into
    the driver's multiprocessing workers.
    """
    fn_info = src_tmpl % ('info', 'csv')
    reduce_info = pd.read_csv(fn_info).iloc[0].to_dict()
    sim_dir = os.path.dirname(fn_info) or '.'
    reduce_info.update(PlanetNames.from_dir(sim_dir).to_reduce_info())
    return reduce_info
