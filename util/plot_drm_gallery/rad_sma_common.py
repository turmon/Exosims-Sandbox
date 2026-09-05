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
from scipy.stats import gaussian_kde
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

## Kernel-density estimates over this plane (see kde_on_bins)
# Grid points per axis.  The estimate costs O(n_points x grid^2) and is all of
# the run time of the plots that use it, so the grid is worth sizing: the
# kernel bandwidth comes out around 0.1 dex against a plot spanning ~2.5 x 1.5
# dex, so 96 points still samples one bandwidth about 5 times.  Measured
# against grid=160, that changes the density by under half a percent, for a
# third of the time.
KDE_GRID = 96
# Sample cap.  The planet-population table runs to ~400k rows for a 100-DRM
# ensemble, which would take minutes; 10k points is already far more than a
# 2-D density needs.  Unlike the grid, this one is not free to reduce --
# halving it moves the surface by ~12% of peak, which is Monte-Carlo noise.
KDE_MAX_POINTS = 10000
# contour levels, as a fraction of the peak density: below the first one,
# nothing is filled, so whatever is drawn underneath stays visible
KDE_LEVELS = np.linspace(0.05, 1.0, 10)
# Ratio maps (kde_ratios_on_bins) are probabilities, so their levels are
# absolute, and the same for every such plot: the colors mean one thing.
RATIO_LEVELS = np.linspace(0.0, 1.0, 11)
# Where the denominator has essentially no data the ratio is meaningless, so
# it is masked: cells holding less than this fraction of the peak kernel
# weight come back NaN, and contourf leaves them blank.
RATIO_WEIGHT_FLOOR = 1e-3
# grid points evaluated per chunk, to bound the (n_points x chunk) work array
RATIO_CHUNK = 512


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


def kde_on_bins(sma, rp, binner, grid=KDE_GRID, max_points=KDE_MAX_POINTS, seed=0):
    r'''Kernel density of points in the radius/SMA plane, over the 5x3 bins.

    Estimated in log10 coordinates, because the plane is plotted log-log: a
    Gaussian kernel in linear SMA would be badly mis-shaped at the low end.
    The returned grid is in data coordinates, ready to hand to contourf, and
    the density is per dex^2.

    Large samples are randomly subsampled to max_points -- with a fixed seed,
    so re-running reproduces the plot -- because the cost is the product of
    the sample size and the grid size.  Returns (X, Y, Z, n_used).

    Raises ValueError or numpy.linalg.LinAlgError if the sample is too small
    or degenerate for a density estimate; the caller decides what to say.
    '''
    x, y = np.log10(sma), np.log10(rp)
    n_used = len(x)
    if n_used > max_points:
        inx = np.random.default_rng(seed).choice(n_used, max_points, replace=False)
        x, y = x[inx], y[inx]
        n_used = max_points
    kernel = gaussian_kde(np.vstack((x, y)))
    (sma_lo, sma_hi), (rp_lo, rp_hi) = koppa_bin_extent(binner)
    xg = np.linspace(np.log10(sma_lo), np.log10(sma_hi), grid)
    yg = np.linspace(np.log10(rp_lo),  np.log10(rp_hi),  grid)
    Xg, Yg = np.meshgrid(xg, yg)
    Z = kernel(np.vstack((Xg.ravel(), Yg.ravel()))).reshape(Xg.shape)
    return 10.0**Xg, 10.0**Yg, Z, n_used


def kde_ratios_on_bins(sma, rp, masks, binner, grid=KDE_GRID,
                           max_points=KDE_MAX_POINTS, seed=0):
    r'''Kernel-regression estimate of P(mask | position) over the 5x3 bins.

    Given points (sma, rp) and one or more boolean masks over those same
    points, return P(mask is true | this position) as a smooth map, for each
    mask, along with the grid to plot it on.

    This is deliberately *not* a ratio of two separately-fitted densities.
    Two gaussian_kde fits choose their bandwidth and orientation from their
    own samples, so their ratio is not bounded, blows up wherever the
    denominator thins out, and is not a probability.  Instead the numerator
    and denominator here are kernel sums over the *same* sample with the
    *same* kernel:

        P(mask | x) = sum_i w_i(x) mask_i / sum_i w_i(x)

    which is the Nadaraya-Watson estimator of the indicator, and lies in
    [0, 1] by construction.  Cells where the denominator holds less than
    RATIO_WEIGHT_FLOOR of its peak weight come back NaN rather than as a
    ratio of two nearly-zero numbers.

    The kernel is the one gaussian_kde would choose for the denominator
    sample, so these maps and the density maps are smoothed alike.  Several
    masks are evaluated in one pass because they share that kernel.

    Returns (X, Y, {name: R}, n_used).
    '''
    x, y = np.log10(sma), np.log10(rp)
    n_used = len(x)
    if n_used > max_points:
        inx = np.random.default_rng(seed).choice(n_used, max_points, replace=False)
        x, y = x[inx], y[inx]
        masks = {name: m[inx] for name, m in masks.items()}
        n_used = max_points
    # borrow scipy's bandwidth choice, then do the weighting ourselves
    pts = np.vstack((x, y))
    inv_cov = gaussian_kde(pts).inv_cov

    (sma_lo, sma_hi), (rp_lo, rp_hi) = koppa_bin_extent(binner)
    xg = np.linspace(np.log10(sma_lo), np.log10(sma_hi), grid)
    yg = np.linspace(np.log10(rp_lo),  np.log10(rp_hi),  grid)
    Xg, Yg = np.meshgrid(xg, yg)
    flat = np.vstack((Xg.ravel(), Yg.ravel()))
    n_grid = flat.shape[1]

    den = np.zeros(n_grid)
    num = {name: np.zeros(n_grid) for name in masks}
    weights = {name: m.astype(float) for name, m in masks.items()}
    for lo in range(0, n_grid, RATIO_CHUNK):
        hi = min(lo + RATIO_CHUNK, n_grid)
        # squared Mahalanobis distance from each grid point to each sample
        dx = flat[0, lo:hi, None] - pts[0][None, :]
        dy = flat[1, lo:hi, None] - pts[1][None, :]
        d2 = (inv_cov[0, 0] * dx * dx + 2 * inv_cov[0, 1] * dx * dy
                  + inv_cov[1, 1] * dy * dy)
        w = np.exp(-0.5 * d2)
        den[lo:hi] = w.sum(axis=1)
        for name, wt in weights.items():
            num[name][lo:hi] = w @ wt

    ok = den > den.max() * RATIO_WEIGHT_FLOOR
    ratios = {}
    for name in masks:
        r = np.full(n_grid, np.nan)
        r[ok] = num[name][ok] / den[ok]
        ratios[name] = r.reshape(Xg.shape)
    return 10.0**Xg, 10.0**Yg, ratios, n_used


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
