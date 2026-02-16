#!/usr/bin/env python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
import os
# uniform plot appearance
from . import common_style as cs

# Program name for error messages
PROGNAME = os.path.basename(sys.argv[0])

# Verbosity (also set from mode)
VERBOSE = 1


def plot_drm_fuel_used(src_tmpl, dest_tmpl, mode):
    """
    Plot fuel use versus time in a drm-set
    
    Time-series plots of fuel use. The "t_fuel" input is the same as
    "t_det_time" elsewhere -- it contains more than just fuel.
    
    Parameters
    ----------
    src_tmpl : str
        Template string for input file paths with two %s placeholders
        Example: "data/%s.%s" will be filled with ("info", "csv") and ("fuel", "csv")
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
    
    # Load data using the source template
    # update global verbosity
    global VERBOSE
    VERBOSE = mode.get('verbose', VERBOSE)

    try:
        info_file = src_tmpl % ("info", "csv")
        t_info = pd.read_csv(info_file)
    except Exception as e:
        print(f"Warning: Could not load info file: {e}")
        t_info = pd.DataFrame()
    
    try:
        fuel_file = src_tmpl % ("times", "csv")
        t_fuel = pd.read_csv(fuel_file)
    except Exception as e:
        print(f"Error: Could not load fuel file: {e}")
        sys.exit(1)
    
    # File extensions to write
    ext_list = ['png']
    
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
        return
    
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
    title2 = cs.plot_make_title(t_info)
    
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
    if dest_tmpl:
        for ext in ext_list:
            fn_gfx = dest_tmpl % ('fuel', ext)
            if VERBOSE:
                print(f'\tExport: {fn_gfx}')
            # figure background is transparent, axes not transparent
            fig.patch.set_facecolor('none')
            fig.savefig(fn_gfx, dpi=200, bbox_inches='tight')
    
    plt.close(fig)
    
    ####################################################################
    # Delta-v plot
    ####################################################################
    
    # Check if delta-v data exists
    if 'h_time_delta_v_slew_cume_mean' not in t_fuel.columns:
        print('No delta-v in the given table, skipping. Redo reduce to fix.')
        return
    
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
    title2 = cs.plot_make_title(t_info)
    
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
    if dest_tmpl:
        for ext in ext_list:
            fn_gfx = dest_tmpl % ('delta-v', ext)
            if VERBOSE:
                print(f'\tExport: {fn_gfx}')
            # figure background is transparent, axes not transparent
            fig.patch.set_facecolor('none')
            fig.savefig(fn_gfx, dpi=200, bbox_inches='tight')
    
    plt.close(fig)


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
    parser.add_argument('--mode_op', type=str, default='*',
                       help='Operation mode (default: "*")')
    
    args = parser.parse_args()
    
    # Create mode dictionary
    mode = {'op': args.mode_op}
    
    # Run the plotting function
    plot_drm_fuel_used(args.src_tmpl, args.dest_tmpl, mode)
    

if __name__ == '__main__':
    main()
