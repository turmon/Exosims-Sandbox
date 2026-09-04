#!/usr/bin/env python
r"""
Plot successful characterizations in radius/SMA ("Kopparapu") coordinates.

One point per successful characterization of an earthlike planet, pooled over
every simulation in the ensemble, on the same radius vs. luminosity-scaled-SMA
axes as the binned rate tables.  Drawn twice: as a scatter, and as a kernel
density estimate.

The 5x3 radius/insolation bins and the earthlike region are drawn behind the
data, so a point can be located against the class definition it satisfies.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
import os
# uniform plot appearance
try:
    from . import common_style as cs
    from . import rad_sma_common as rsc
except ImportError:
    import common_style as cs
    import rad_sma_common as rsc

# Program name for error messages
PROGNAME = os.path.basename(sys.argv[0])

# columns we cannot do without (added to reduce-earth-char-list.csv in 9/2026)
NEEDED_COLUMNS = ('is_success', 'Rp', 'sma_scaled')


def plot_drm_rad_sma_chars(reduce_info, plot_data, dest_tmpl, mode):
    """
    Plot successful characterizations in radius/SMA coordinates.

    Parameters
    ----------
    reduce_info : dict
        Metadata dict from reduce-info.csv
    plot_data : list of DataFrame
        Pre-loaded CSV data [earth-char-list]
    dest_tmpl : str
        Template string for output file paths with two %s placeholders
    mode : dict
        Dictionary with 'op' key containing operation mode string

    Outputs
    -------
    Saves plots to disk with names:
        rad-sma-chars-scatter.png
        rad-sma-chars-density.png
    """

    # Unpack CSV data
    t_earth_chars, = plot_data

    # display names of the earthlike planet class (may be re-defined
    # in config-reduce.json, e.g. to a Sub-Neptune population)
    pn = cs.planet_names(reduce_info)

    # Missions without chars have a near-empty table; older reductions have no
    # Rp/sma_scaled columns.  Neither is an error.
    missing = [c for c in NEEDED_COLUMNS if c not in t_earth_chars.columns]
    if missing:
        print(f'\t{PROGNAME}: No {", ".join(missing)} in the {pn.short} char list, '
                  'skipping (re-run reduction?)')
        return []

    # Successful chars only, pooled over the whole ensemble.
    # Coerce rather than trust the dtype: a table with no rows reads back as
    # object-dtype columns, on which the tests below would raise rather than
    # come out empty.  Anything non-numeric becomes NaN and drops out.
    def column(name):
        return pd.to_numeric(t_earth_chars[name], errors='coerce').values

    ok = column('is_success') > 0
    sma = column('sma_scaled')[ok]
    rp = column('Rp')[ok]
    # guard against non-positive values: these axes are log-log
    good = (sma > 0) & (rp > 0) & np.isfinite(sma) & np.isfinite(rp)
    sma, rp = sma[good], rp[good]
    n_ok = len(sma)
    if n_ok == 0:
        print(f'\t{PROGNAME}: No successful {pn.short} chars, skipping')
        return []

    # bin geometry, as customized by config-reduce.json for this scenario
    binner = rsc.configured_binner(reduce_info.get('_sim_dir', '.'), log_origin=PROGNAME)

    # Track output files
    tracker = cs.PlotTracker(ext_list=mode.get('ext_list'), reduce_info=reduce_info)

    # Inner function: Set up plot/axis styles, title, axis labels
    def style_rad_sma_plot(ax, title1):
        """Style the radius/SMA plot with title, labels, and limits"""
        # Format the title
        title2 = cs.plot_make_title(reduce_info)

        # Set title (preventing special interpretation of _) with bold
        ax.set_title(f'{title2}\n{title1}', fontsize=11*1.1, fontweight='bold')

        # Axis labels
        ax.set_xlabel('Semi-Major Axis, Luminosity-Scaled [AU]', fontweight='bold')
        ax.set_ylabel('$R_p$ [Earth radii]', fontweight='bold')
        ax.tick_params(labelsize=13)
        # frame the whole 5x3 grid, log-log, and un-exponentiate the tick labels
        rsc.set_koppa_limits(ax, binner)
        for axis in (ax.xaxis, ax.yaxis):
            axis.set_major_formatter(plt.FuncFormatter(lambda v, _: '{:.8g}'.format(v)))
        ax.grid(True, alpha=0.3)

    # Inner function: write the current figure to files
    def write_plots(fig, dest_name):
        """Write the current figure to various files"""
        tracker.write_plots(fig, dest_name, dest_tmpl, verbose=mode.get('verbose', 1))

    ####################################################################
    # A: Scatter
    ####################################################################

    fig, ax = plt.subplots(figsize=(8.5, 5))
    rsc.draw_koppa_boxes(ax, binner)
    rsc.draw_earthlike_region(ax, binner)
    # small, translucent: these pile up on top of each other
    ms = 4 if n_ok < 500 else 1
    ax.plot(sma, rp, linestyle='none', marker='.', markersize=ms,
            color='black', alpha=0.35, zorder=3,
            label=f'Successful chars (N = {n_ok})')
    style_rad_sma_plot(
        ax, f'Successful {pn.name} Characterizations, Pooled Over Ensemble')
    ax.legend(loc='upper left', framealpha=0.8)
    write_plots(fig, 'rad-sma-chars-scatter')
    plt.close(fig)

    ####################################################################
    # B: Kernel density
    ####################################################################

    # estimated in log coordinates, on the bin grid -- see rsc.kde_on_bins
    try:
        Xg, Yg, Z, _n_used = rsc.kde_on_bins(sma, rp, binner)
    except (ValueError, np.linalg.LinAlgError) as e:
        # too few points, or all of them collinear/identical
        # (scipy's message is a paragraph; the first sentence is the reason)
        print(f'\t{PROGNAME}: Skipping {pn.short} char density plot '
                  f'({n_ok} points): {str(e).split(".")[0]}')
        return tracker.get_files()

    fig, ax = plt.subplots(figsize=(8.5, 5))
    # bins are context here, not the subject: hold them back
    rsc.draw_koppa_boxes(ax, binner, alpha=0.30)
    cs_kde = ax.contourf(Xg, Yg, Z, levels=rsc.KDE_LEVELS * Z.max(),
                         cmap='magma_r', alpha=0.85, zorder=3, extend='max')
    # over the density, not under it: the kernel smooths across the class
    # boundary, so the boundary has to stay visible to read the plot
    rsc.draw_earthlike_region(ax, binner, zorder=4)
    style_rad_sma_plot(
        ax, f'Successful {pn.name} Characterization Density ({n_ok} chars)')
    cbar = fig.colorbar(cs_kde, ax=ax)
    cbar.set_label('Probability Density [/ dex$^2$]', fontweight='bold')
    write_plots(fig, 'rad-sma-chars-density')
    plt.close(fig)

    return tracker.get_files()


def main():
    r"""
    Command-line interface for plot_drm_rad_sma_chars
    """
    parser = argparse.ArgumentParser(
        description='Plot successful characterizations in radius/SMA coordinates',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
    plot_drm_rad_sma_chars.py SOURCE_TEMPLATE OUTPUT_TEMPLATE

The first argument is the source template with two %s placeholders that will be
filled with ("info", "csv") and ("earth-char-list", "csv") to locate input files.

The second argument is the destination template with two %s placeholders for
the plot name and file extension.
        """
    )

    parser.add_argument('src_tmpl', type=str,
                       help='Source template string (e.g., "sims/scenario/reduce-%%s.%%s")')
    parser.add_argument('dest_tmpl', type=str,
                       help='Destination template string (e.g., "sims/scenario/gfx/det-%%s.%%s")')
    parser.add_argument('--mode_op', type=str, default='',
                       help='Operation mode, default: "" (normal)')
    parser.add_argument('--verbose', '-v', action='count', default=1,
                       help='Verbosity')
    parser.add_argument('--quiet', '-q', action='store_true', help='Minimal verbosity')

    args = parser.parse_args()
    if args.quiet: args.verbose = 0

    # Create mode dictionary
    mode = {'op': args.mode_op, 'verbose': args.verbose}

    # Read info file and convert to dict (plus planet-class display names)
    reduce_info = cs.load_reduce_info(args.src_tmpl)

    # Load CSV data and run the plotting function
    plot_data = cs.load_csv_files(args.src_tmpl, ['earth-char-list'])
    rv = plot_drm_rad_sma_chars(reduce_info, plot_data, args.dest_tmpl, mode)
    return rv


if __name__ == '__main__':
    rv = main()
    if rv is None:
        print(f"Plots failed. Error signaled.", file=sys.stderr)
    else:
        print(f"Done. Wrote {len(rv)} plot(s).")
    sys.exit(1 if rv is None else 0)
