#!/usr/bin/env python
"""
Common styles, functions, etc.

"""

import pandas as pd


class PlotTracker:
    """Track graphics files written by a plot routine."""

    def __init__(self):
        self._files = []               # list of filenames written
        self._ext_list = ['png']       # default

    def set_ext_list(self, ext_list):
        """Set the file extensions for subsequent writes."""
        self._ext_list = list(ext_list)

    def write_plots(self, fig, dest_name, dest_tmpl,
                    ext_list=None, verbose=True, dpi=200, facecolor='none'):
        """Write figure to file(s) and record what was written."""
        if ext_list is None:
            ext_list = self._ext_list
        if dest_tmpl:
            for ext in ext_list:
                fn_gfx = dest_tmpl % (dest_name, ext)
                if verbose:
                    print(f'\tExport: {fn_gfx}')
                if facecolor is not None:
                    fig.patch.set_facecolor(facecolor)
                fig.savefig(fn_gfx, dpi=dpi, bbox_inches='tight')
                self._files.append(fn_gfx)

    def get_files(self):
        """Return list of filenames written."""
        return list(self._files)


# Helper function to create title from t_info
def plot_make_title(t_info):
    """Create plot title from metadata"""
    if t_info.empty:
        rv = ''
    elif 'experiment' in t_info:
        exp_name = str.strip(t_info['experiment'].iloc[0])
        if len(exp_name) < 50:
            chaser = f", Ensemble Size {t_info['ensemble_size'].iloc[0]}"
        else:
            chaser = ''
        rv = exp_name + chaser
    else:
        rv = ''
    return rv
    
