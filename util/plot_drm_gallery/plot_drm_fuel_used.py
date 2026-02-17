#!/usr/bin/env python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
import os
# uniform plot appearance
try:
    from . import common_style as cs
except ImportError:
    import common_style as cs

# Program name for error messages
PROGNAME = os.path.basename(sys.argv[0])


def plot_drm_fuel_used(reduce_info, plot_data, dest_tmpl, mode):
    """
    Plot fuel use versus time in a drm-set

    Time-series plots of fuel use. The "t_fuel" input is the same as
    "t_det_time" elsewhere -- it contains more than just fuel.

    Parameters
    ----------
    reduce_info : dict
        Metadata dict from reduce-info.csv
    plot_data : list of DataFrame
        Pre-loaded CSV data [times]
    dest_tmpl : str
        Template string for output file paths (should contain {} for name and extension)
    mode : dict
        Dictionary with operation mode settings

    Outputs
    -------
    Saves plots to disk with names:
        fuel.png
        delta-v.png
    """

    # Unpack CSV data
    t_fuel, = plot_data

    # Track output files
    tracker = cs.PlotTracker(ext_list=mode.get('ext_list'))

    # Inner function: write the current figure to files
    def write_plots(fig, dest_name):
        """Write the current figure to various files"""
        tracker.write_plots(fig, dest_name, dest_tmpl, verbose=mode['verbose'])

    ####################################################################
    # Setup
    ####################################################################
    
    # Manual line color order
    # purple is for slew, black for total
    line_colors = ['tab:purple', 'tab:orange', 'black']

    # uniform error bar properties
    ebar_props = {"marker": '.',
                  "linewidth": 1.7,
                  "errorevery": 4,
                  "elinewidth": 1,
                  "capsize": 1.5}


    ####################################################################
    # Fuel (in kg) plot
    ####################################################################
    
    # Check if fuel data exists
    if 'h_time_fuel_slew_mean' not in t_fuel.columns:
        print('No fuel use in the given table, skipping.')
        return []
    
    # Create figure
    fig, ax = plt.subplots(figsize=(8.5, 5))
    
    # What to plot
    names = ['h_time_fuel_slew', 'h_time_fuel_keep', 'h_time_fuel_all']
    names_legend = ['Slew', 'Stationkeeping', 'Total']
    t_offsets = [0, 10, 5]
    n_plot = len(names)
    
    # Sample times
    tsamp = t_fuel['h_det_time_lo'].values
    
    # Put each of the above on one plot
    top = 0
    for n, name in enumerate(names):
        f_mean = f'{name}_mean'
        f_std = f'{name}_std'
        
        ax.errorbar(tsamp + t_offsets[n],
                    t_fuel[f_mean].values,
                    yerr=t_fuel[f_std].values,
                    color=line_colors[n],
                    **ebar_props)
        
        # Get some control over axis range
        top = max(top, t_fuel[f_mean].max())
    
    # Plot/axis styles
    title1 = 'Cumulative Fuel Use vs. Time'
    title2 = cs.plot_make_title(reduce_info)
    
    # Control axis y-range:
    # negative not allowed
    # top of range can't be too big
    # (if top guards against zero fuel use for empty DRMs)
    if top > 0:
        ylim = ax.get_ylim()
        ax.set_ylim(max(0, ylim[0]), 1.2 * top)
    
    # Title: prevent special interpretation of _ and make bold
    ax.set_title(f'{title2}\n{title1}', fontsize=11*1.1, fontweight='bold')
    ax.set_xlabel('Time [days]', fontweight='bold')
    ax.set_ylabel('Fuel Used [kg]', fontweight='bold')
    ax.legend(names_legend, loc='upper left')
    ax.tick_params(labelsize=13)
    ax.grid(True)
    
    # Write the plot out
    write_plots(fig, 'fuel')

    plt.close(fig)
    
    ####################################################################
    # Delta-v plot
    ####################################################################
    
    # Check if delta-v data exists
    if 'h_time_delta_v_slew_cume_mean' not in t_fuel.columns:
        print('No delta-v in the given table, skipping. Redo reduce to fix.')
        return tracker.get_files()
    
    # Create figure
    fig, ax = plt.subplots(figsize=(8.5, 5))
    
    # What to plot
    # "total" makes no sense for starshades (two separate objects moving),
    # so do not display the sum (.._all_...) as a separate line.
    names = ['h_time_delta_v_slew_cume', 'h_time_delta_v_obs_cume']
    names_legend = ['Slew', 'Observing']
    t_offsets = [0, 10, 5]
    n_plot = len(names)
    
    # Sample times
    tsamp = t_fuel['h_det_time_lo'].values
    
    # Put each of the above on one plot
    top = 0
    for n, name in enumerate(names):
        f_mean = f'{name}_mean'
        f_std = f'{name}_std'
        
        ax.errorbar(tsamp + t_offsets[n],
                    t_fuel[f_mean].values,
                    yerr=t_fuel[f_std].values,
                    color=line_colors[n],
                    **ebar_props)
        
        # Get some control over axis range
        top = max(top, t_fuel[f_mean].max())
    
    # Plot/axis styles
    title1 = 'Cumulative Delta-V vs. Time'
    title2 = cs.plot_make_title(reduce_info)
    
    # Control axis y-range:
    # negative not allowed
    # top of range can't be too big
    # (if top guards against zero fuel use for empty DRMs)
    if top > 0:
        ylim = ax.get_ylim()
        ax.set_ylim(max(0, ylim[0]), 1.2 * top)
    
    # Title: prevent special interpretation of _ and make bold
    ax.set_title(f'{title2}\n{title1}', fontsize=11*1.1, fontweight='bold')
    ax.set_xlabel('Time [days]', fontweight='bold')
    ax.set_ylabel('Delta V [m/s]', fontweight='bold')
    ax.legend(names_legend, loc='upper left')
    ax.tick_params(labelsize=13)
    ax.grid(True)
    
    # Write the plot out
    write_plots(fig, 'delta-v')

    plt.close(fig)

    return tracker.get_files()


def main():
    r"""
    Command-line interface for plot_drm_fuel_used
    """
    parser = argparse.ArgumentParser(
        description='Plot fuel use versus time in a DRM-set',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
    plot_drm_fuel_used.py SOURCE_TEMPLATE OUTPUT_TEMPLATE
    
The first argument is the source template with two %s placeholders that will be
filled with ("info", "csv") and ("time", "csv") to locate input files.

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

    # Read info file and convert to dict
    info_file = args.src_tmpl % ("info", "csv")
    reduce_info = pd.read_csv(info_file).iloc[0].to_dict()

    # Load CSV data and run the plotting function
    plot_data = cs.load_csv_files(args.src_tmpl, ['times'])
    rv = plot_drm_fuel_used(reduce_info, plot_data, args.dest_tmpl, mode)
    return rv


if __name__ == '__main__':
    rv = main()
    if rv is None:
        print(f"Plots failed. Error signaled.", file=sys.stderr)
    else:
        print(f"Done. Wrote {len(rv)} plot(s).")
    sys.exit(1 if rv is None else 0)
