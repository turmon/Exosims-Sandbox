#!/usr/bin/env python
"""
Common styles, functions, etc.

"""

import os
import sys
import pandas as pd


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
                print(f'\tExport: {fn_gfx}')
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

