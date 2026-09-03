#!/usr/bin/env python
"""
Common styles, functions, etc.

"""

import os
import sys
import pandas as pd
import matplotlib as mpl

# Program name for error messages
PROGNAME = os.path.basename(sys.argv[0])

# reduce_drm_tools lives one level up, in util/.  Make it importable whether we
# were started by the driver (which puts util/ on sys.path) or standalone from
# within this directory.
_UTIL_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _UTIL_DIR not in sys.path:
    sys.path.append(_UTIL_DIR)
from reduce_drm_tools.PlanetNames import PlanetNames

# Ensemble-size annotation: which end of the plot box it sits at ('left' or
# 'right'), and the fallback corner, in figure coordinates
ENSEMBLE_NOTE_SIDE = 'right'
ENSEMBLE_NOTE_MARGIN = 0.01
ENSEMBLE_NOTE_BOTTOM = 0.00


def _ensemble_note_corner():
    """Fallback position and alignment: bottom corner on the chosen side."""
    if ENSEMBLE_NOTE_SIDE == 'right':
        return 1.0 - ENSEMBLE_NOTE_MARGIN, ENSEMBLE_NOTE_BOTTOM, 'right'
    return ENSEMBLE_NOTE_MARGIN, ENSEMBLE_NOTE_BOTTOM, 'left'


def plot_ensemble_note(reduce_info):
    """Return the ensemble-size annotation text ('' if the size is unknown)."""
    if not reduce_info or 'ensemble_size' not in reduce_info:
        return ''
    return f"(N = {reduce_info['ensemble_size']} runs)"


def _place_ensemble_note(fig, txt):
    """Move the note to the xlabel's baseline, aligned with the plot box edge.

    We write with bbox_inches='tight', which crops to the union of all the
    artists.  A note in a figure corner sits outside that union, so the crop
    grows to reach it and the plot is padded with an empty strip -- up to 7%
    of the width on the wider (9.5in) figures.

    Two moves avoid it.  Horizontally, align to one of the axes' spines,
    which are inboard of the axis labels, tick labels and any colorbar, so
    the note can never be the outermost artist.  Vertically, sit on the
    xlabel's own baseline: the xlabel is centered, so the space beside it is
    already inside the crop.  Net cost in both directions is then about zero.

    Needs a draw to find out where the xlabel landed.  If the xlabel is long
    enough that the note would run into it, we keep the spine alignment
    (which is what saves the width) and drop the note to the bottom of the
    figure instead.  Returns False only if we could not place it at all.
    """
    axes = fig.get_axes()
    if not axes:
        return False
    ax = axes[0]
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    label = ax.xaxis.get_label()
    bb_label = label.get_window_extent(renderer)
    on_right = ENSEMBLE_NOTE_SIDE == 'right'
    box = ax.get_position()
    x = box.x1 if on_right else box.x0
    y = fig.transFigure.inverted().transform((0, bb_label.y0))[1]
    txt.set_position((x, y))
    # keep clear of the xlabel: about half a character of space
    pad = 0.5 * txt.get_fontsize() * fig.dpi / 72.0
    bb_note = txt.get_window_extent(renderer)
    clear = bb_note.x0 > bb_label.x1 + pad if on_right else bb_note.x1 < bb_label.x0 - pad
    if label.get_text() and not clear:
        # no room on the xlabel's line -- go below it, still spine-aligned
        txt.set_position((x, ENSEMBLE_NOTE_BOTTOM))
    return True


def plot_add_ensemble_note(fig, reduce_info):
    """Stamp the ensemble size along the bottom of the figure.

    This information used to be a ", Ensemble Size N" chaser on the title
    line (see plot_make_title).  Which end it goes to is set by
    ENSEMBLE_NOTE_SIDE.  It is placed in *figure* coordinates, not axes
    coordinates, so it lands in the same spot regardless of how the axes
    within the figure are laid out.

    Styled to match the axis labels, but not bold -- it is a footnote, not
    part of the plot proper.

    Idempotent: a figure that is written more than once (under different
    names) is stamped only once.
    """
    text = plot_ensemble_note(reduce_info)
    if not text or getattr(fig, '_has_ensemble_note', False):
        return
    x_corner, y_corner, halign = _ensemble_note_corner()
    txt = fig.text(x_corner, y_corner, text,
                   ha=halign, va='bottom',
                   fontsize=mpl.rcParams['axes.labelsize'],
                   fontweight='normal')
    # tuck it in beside the xlabel if we can -- see _place_ensemble_note
    try:
        if not _place_ensemble_note(fig, txt):
            txt.set_position((x_corner, y_corner))
    except Exception:
        # any backend that will not give us a renderer: corner it is
        txt.set_position((x_corner, y_corner))
    fig._has_ensemble_note = True


class PlotTracker:
    """Track graphics files written by a plot routine."""

    def __init__(self, ext_list=None, reduce_info=None):
        # list of filenames written
        self._files = []               
        # metadata for the lower-left ensemble-size annotation (may be None)
        self._reduce_info = reduce_info
        # file extensions to make
        if ext_list is None:
            self._ext_list = ['png'] # default
        else:
            self._ext_list = ext_list[:]

    def set_ext_list(self, ext_list):
        """Set the file extensions for subsequent writes."""
        self._ext_list = list(ext_list)

    def write_plots(self, fig, dest_name, dest_tmpl,
                    ext_list=None, verbose=1, dpi=200, facecolor='none'):
        """Write figure to file(s) and record what was written.

        verbose is a level, not a flag: 0 is silent, 1 and up name each file
        written.  It is the same scale as mode['verbose'] in the plot
        modules, which is where callers get it.
        """
        if ext_list is None:
            ext_list = self._ext_list
        # uniform ensemble-size annotation on every plot we write
        plot_add_ensemble_note(fig, self._reduce_info)
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
        # the ensemble size used to be appended here; it is now annotated in
        # the plot's lower-left corner (see plot_add_ensemble_note)
        rv = str.strip(reduce_info['experiment'])
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
            print(f"{PROGNAME}: Fatal: Could not load CSV file '{csv_path}': {e}",
                  file=sys.stderr)
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
