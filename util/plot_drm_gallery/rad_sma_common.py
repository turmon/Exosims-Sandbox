#!/usr/bin/env python
r"""
rad_sma_common.py -- shared furniture for radius/SMA ("Kopparapu") plots

The planet radius vs. semi-major axis plane is drawn by two routines that do
otherwise unrelated things:

  util/rad-sma-rectangle-bin-plot.py   -- the binned table of rates/counts
  plot_drm_rad_sma_chars.py            -- characterization attempts, as points

Both need the same 5x3 bin rectangles, the same colors, and the same outline
of the earthlike region, so those live here and not in either one.

Note the x coordinate throughout is the *luminosity-scaled* SMA, a/sqrt(L),
which is what the bins and RpLBins.is_earthlike() are defined in terms of
(SMA = 1/sqrt(L)).  Raw orbital SMA does not belong on these axes.

turmon sep 2026
"""

import os
import sys
import math
import numpy as np
from pathlib import Path
from matplotlib.patches import Rectangle, Polygon

# reduce_drm_tools lives one level up, in util/ -- same dance as common_style
_UTIL_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _UTIL_DIR not in sys.path:
    sys.path.append(_UTIL_DIR)
from reduce_drm_tools.PlanetBins import RpLBins
from reduce_drm_tools import utils


## Style of the polygon outlining the earthlike region
EARTH_STYLE = dict(fill=False,
                   facecolor=None,
                   edgecolor='lightgreen',
                   linewidth=2.0,
                   alpha=0.7, # tiny bit of transparency
                   hatch='/')

## Fill colors of the three luminosity bins, in hot -> cold order
L_BIN_COLORS = ('xkcd:pastel red', 'dodgerblue', 'lightskyblue')


def configured_binner(sim_dir, log_origin=None):
    r'''Return an RpLBins reflecting the config-reduce.json reachable from sim_dir.

    RpLBins refuses to instantiate until customize_parameters() has been
    called, because class-level customization does not survive being forked
    or pickled into a worker process.  So each process that wants a binner
    re-applies it; that is what this does.  An absent config is normal, and
    leaves the default bins.'''
    config = utils.load_reduce_config(Path(sim_dir), log_origin=log_origin) or {}
    RpLBins.customize_parameters(config)
    return RpLBins()


def earthlike_polygon_xy(binner):
    r'''Return the vertices of the "Nevada-shaped" earthlike region, as (x, y).'''
    x = np.array([binner.Earth_SMA_lo, binner.Earth_SMA_hi,
                  binner.Earth_SMA_hi, binner.Earth_SMA_lo])
    y = np.array([binner.Earth_Rp_hi,  binner.Earth_Rp_hi,
                  binner.Earth_Rp_lo2, binner.Earth_Rp_lo1])
    return x, y


def koppa_bin_extent(binner):
    r'''Return ((sma_lo, sma_hi), (Rp_lo, Rp_hi)) spanned by the 5x3 bins.'''
    sma = 1.0 / np.sqrt(np.asarray(binner.L_bins, dtype=float))
    return (sma.min(), sma.max()), (binner.Rp_bins[0], binner.Rp_bins[-1])


def draw_koppa_boxes(ax, binner, alpha=1.0, zorder=0):
    r'''Draw the 5x3 radius/insolation bins as colored rectangles.

    Rp_bins gives the y (radius) stripes; within each stripe, L_bins gives the
    x (SMA) cuts, converted by SMA = 1/sqrt(L).  Pass a reduced alpha when the
    boxes are background for something drawn on top of them.'''
    for r_inx in range(len(binner.Rp_bins) - 1):
        r0, r1 = binner.Rp_bins[r_inx:r_inx+2]
        # L bins run hot -> cold, so the SMAs come out in increasing order
        sma_bins = [1.0 / math.sqrt(L) for L in binner.L_bins[r_inx]]
        for L_inx in range(len(sma_bins) - 1):
            a0, a1 = sma_bins[L_inx:L_inx+2]
            ax.add_patch(Rectangle((a0, r0), a1 - a0, r1 - r0,
                                   facecolor=L_BIN_COLORS[L_inx],
                                   edgecolor='white',
                                   alpha=alpha, zorder=zorder))


def draw_earthlike_region(ax, binner, zorder=2):
    r'''Outline the earthlike region, in the style used across the Sandbox.'''
    x, y = earthlike_polygon_xy(binner)
    ax.add_patch(Polygon(np.vstack((x, y)).T, zorder=zorder, **EARTH_STYLE))


def set_koppa_limits(ax, binner, margin=0.03):
    r'''Set log-log limits framing the 5x3 bins, with a small margin in dex.'''
    (x0, x1), (y0, y1) = koppa_bin_extent(binner)
    def _span(lo, hi):
        pad = margin * (math.log10(hi) - math.log10(lo))
        return 10.0**(math.log10(lo) - pad), 10.0**(math.log10(hi) + pad)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(*_span(x0, x1))
    ax.set_ylim(*_span(y0, y1))
