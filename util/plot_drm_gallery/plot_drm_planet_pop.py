#!/usr/bin/env python
r"""
Plot the simulated planet population, and its yield, as densities.

Three kernel-density plots over the radius vs. luminosity-scaled-SMA plane,
from reduce-planet-population.csv: all the planets in the table, those that
were detected, and those that were characterized.  Read together they are the
mission's funnel, in the plane where the planet classes are defined.

The table holds planets around stars the mission visited, so "all" is the
population it had the chance to observe, not the whole simulated universe.
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

# columns we cannot do without
NEEDED_COLUMNS = ('sma_scaled', 'radius', 'det_ok', 'char_ok')

# the three subsets: (file stem, selector column or None, title phrase)
SUBSETS = (
    ('all',  None,      'All Planets'),
    ('det',  'det_ok',  'Detected Planets'),
    ('char', 'char_ok', 'Characterized Planets'),
    )


def plot_drm_planet_pop(reduce_info, plot_data, dest_tmpl, mode):
    """
    Plot planet-population densities in radius/SMA coordinates.

    Parameters
    ----------
    reduce_info : dict
        Metadata dict from reduce-info.csv
    plot_data : list of DataFrame
        Pre-loaded CSV data [planet-population]
    dest_tmpl : str
        Template string for output file paths with two %s placeholders
    mode : dict
        Dictionary with 'op' key containing operation mode string

    Outputs
    -------
    Saves plots to disk with names:
        planet-pop-density-all.png
        planet-pop-density-det.png
        planet-pop-density-char.png
    """

    # Unpack CSV data
    t_planets, = plot_data

    # display names of the earthlike planet class (may be re-defined
    # in config-reduce.json, e.g. to a Sub-Neptune population)
    pn = cs.planet_names(reduce_info)

    # An older reduction has no planet-population table at all, in which case
    # the driver skips us before this; guard the column set regardless.
    missing = [c for c in NEEDED_COLUMNS if c not in t_planets.columns]
    if missing:
        print(f'\t{PROGNAME}: No {", ".join(missing)} in the planet population, '
                  'skipping (re-run reduction?)')
        return []

    # Coerce rather than trust the dtype: a table with no rows reads back as
    # object-dtype columns.  Anything non-numeric becomes NaN and drops out.
    def column(name):
        return pd.to_numeric(t_planets[name], errors='coerce').values

    sma_all, rp_all = column('sma_scaled'), column('radius')
    # these axes are log-log, so non-positive and non-finite values cannot plot
    good = (sma_all > 0) & (rp_all > 0) & np.isfinite(sma_all) & np.isfinite(rp_all)
    if not np.any(good):
        print(f'\t{PROGNAME}: No usable planet rows, skipping')
        return []

    # bin geometry, as customized by config-reduce.json for this scenario
    binner = rsc.configured_binner(reduce_info.get('_sim_dir', '.'), log_origin=PROGNAME)

    # Track output files
    tracker = cs.PlotTracker(ext_list=mode.get('ext_list'), reduce_info=reduce_info)

    # Inner function: Set up plot/axis styles, title, axis labels
    def style_rad_sma_plot(ax, title1):
        """Style the radius/SMA plot with title, labels, and limits"""
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
    # One density per subset: all planets, detected, characterized
    ####################################################################

    for stem, selector, phrase in SUBSETS:
        take = good if selector is None else (good & (column(selector) > 0))
        sma, rp = sma_all[take], rp_all[take]
        n_pl = len(sma)
        if n_pl == 0:
            print(f'\t{PROGNAME}: No {phrase.lower()} in the population, '
                      f'skipping the {stem} density plot')
            continue
        try:
            Xg, Yg, Z, n_used = rsc.kde_on_bins(sma, rp, binner)
        except (ValueError, np.linalg.LinAlgError) as e:
            # too few points, or all of them collinear/identical
            # (scipy's message is a paragraph; the first sentence is the reason)
            print(f'\t{PROGNAME}: Skipping the {stem} density plot '
                      f'({n_pl} points): {str(e).split(".")[0]}')
            continue

        fig, ax = plt.subplots(figsize=(8.5, 5))
        # bins are context here, not the subject: hold them back
        rsc.draw_koppa_boxes(ax, binner, alpha=0.30)
        cs_kde = ax.contourf(Xg, Yg, Z, levels=rsc.KDE_LEVELS * Z.max(),
                             cmap='magma_r', alpha=0.85, zorder=3, extend='max')
        # over the density, not under it: the kernel smooths across the class
        # boundary, so the boundary has to stay visible to read the plot
        rsc.draw_earthlike_region(ax, binner, zorder=4)
        # say when the density rests on a subsample of the rows
        sampled = '' if n_used == n_pl else f', {n_used} sampled'
        style_rad_sma_plot(ax, f'{phrase}: Density ({n_pl} planets{sampled})')
        cbar = fig.colorbar(cs_kde, ax=ax)
        cbar.set_label('Probability Density [/ dex$^2$]', fontweight='bold')
        write_plots(fig, f'planet-pop-density-{stem}')
        plt.close(fig)

    return tracker.get_files()


def main():
    r"""
    Command-line interface for plot_drm_planet_pop
    """
    parser = argparse.ArgumentParser(
        description='Plot planet-population densities in radius/SMA coordinates',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
    plot_drm_planet_pop.py SOURCE_TEMPLATE OUTPUT_TEMPLATE

The first argument is the source template with two %s placeholders that will be
filled with ("info", "csv") and ("planet-population", "csv") to locate inputs.

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
    plot_data = cs.load_csv_files(args.src_tmpl, ['planet-population'])
    rv = plot_drm_planet_pop(reduce_info, plot_data, args.dest_tmpl, mode)
    return rv


if __name__ == '__main__':
    rv = main()
    if rv is None:
        print(f"Plots failed. Error signaled.", file=sys.stderr)
    else:
        print(f"Done. Wrote {len(rv)} plot(s).")
    sys.exit(1 if rv is None else 0)
