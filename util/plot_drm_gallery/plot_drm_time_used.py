#!/usr/bin/env python
"""
Plot time-used in a drm-set

Time-series plots of time-used by detections, chars, slews, 
both cumulatively over the mission, and month-by-month.

Note: This plot family originally included detection counts (cume,
uniq, revisit) vs. time - within the t_det_time table. We later created 
a newer table (t_yield_time) that has det/char, and allplanets/earth,
and multiband char info. That plot family has superseded some of the
plots that used to be here.
"""

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


def plot_drm_time_used(src_tmpl, dest_tmpl, mode):
    """
    Plot time-used in a drm-set
    
    Time-series plots of time-used by detections, chars, slews, 
    both cumulatively over the mission, and month-by-month.
    
    Parameters
    ----------
    src_tmpl : str
        Template string for input file paths with two %s placeholders
        Example: "data/%s.%s" will be filled with ("info", "csv") and ("det_time", "csv")
    dest_tmpl : str
        Template string for output file paths with two %s placeholders
    mode : dict
        Dictionary with 'op' key containing operation mode string
        
    Outputs
    -------
    Saves plots to disk with names:
        obstime-cume.png
    """
    
    # Load data using the source template
    # update global verbosity
    global VERBOSE
    VERBOSE = mode.get('verbose', VERBOSE)

    try:
        info_file = src_tmpl % ("info", "csv")
        t_info = pd.read_csv(info_file)
    except Exception as e:
        print(f"Warning: Could not load info file: {e}", file=sys.stderr)
        t_info = pd.DataFrame()
    
    try:
        det_time_file = src_tmpl % ("times", "csv")
        t_det_time = pd.read_csv(det_time_file)
    except Exception as e:
        print(f"{PROGNAME}: Fatal: Could not load observation time file: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Skip unless mode.op contains our name or a *
    if '*' not in mode.get('op', '') and 'det_time' not in mode.get('op', ''):
        print('Det time plots: skipping, as directed.')
        return []
    
    # File extensions to write
    ext_list = ['png']

    # Track output files
    tracker = cs.PlotTracker()
    tracker.set_ext_list(ext_list)

    # Inner function: Set up plot/axis styles, title, axis labels
    def style_det_plot(ax, title1, ytext, legtext):
        """Style the detection plot with title, labels, and legend"""
        # Clip Y range at 0
        ylim = ax.get_ylim()
        ax.set_ylim(max(0, ylim[0]), ylim[1])
        
        # Format the title
        title2 = cs.plot_make_title(t_info)
        
        # Set title (preventing special interpretation of _) with bold
        ax.set_title(f'{title2}\n{title1}', fontsize=11*1.1, fontweight='bold')
        
        # Axis labels
        ax.set_xlabel('Time [day]', fontweight='bold')
        ax.set_ylabel(ytext, fontweight='bold')
        
        # Legend
        ax.legend(legtext, loc='upper left')
        ax.tick_params(labelsize=13)
        ax.grid(True)
    
    # Inner function: write the current figure to files
    def write_plots(fig, dest_name):
        """Write the current figure to various files"""
        tracker.write_plots(fig, dest_name, dest_tmpl, verbose=VERBOSE)

    # Offsets on various error-bars in the plot 
    # (units of days with one complete sample per ~30 days)
    t_offsets = [0, 10, 5]
    
    # Sample times
    try:
        tsamp = t_det_time['h_det_time_lo'].values
    except KeyError as e:
        print(f"{PROGNAME}: Fatal: Missing required column in det_time file: {e}", 
              file=sys.stderr)
        sys.exit(1)
    
    # Manual line color order
    # for the obstime plot -- purple for slews
    line_colors = ['tab:blue', 'tab:green', 'tab:purple']
    # line_colors = ['tab:blue', 'tab:green', 'tab:orange']

    # Uniform error bar properties
    ebar_props = {"marker": '.',
                  "linewidth": 1.7,
                  "errorevery": 1,
                  "elinewidth": 1,
                  "capsize": 1.5}
    
    ####################################################################
    # Slews and Characterization vs. time -- Cumulative
    ####################################################################
    
    fig, ax = plt.subplots(figsize=(8.5, 5))
    
    names = ['h_time_det_cume', 'h_time_char_cume', 'h_time_slew_cume']
    names_legend = ['Detection Time', 'Characterization Time', 'Slew Time']
    n_plot = len(names)
    
    # Put each of the above detection-times on one plot
    for n, name in enumerate(names):
        f_mean = f'{name}_mean'
        f_std = f'{name}_std'
        
        ax.errorbar(tsamp + t_offsets[n],
                    t_det_time[f_mean].values,
                    yerr=t_det_time[f_std].values,
                    color=line_colors[n],
                    **ebar_props)
    
    style_det_plot(ax, 
                  'Mission Observation Scheduling (Cumulative) vs. Mission Time',
                  'Cumulative Time Used [day]', 
                  names_legend)
    write_plots(fig, 'obstime-cume')
    plt.close(fig)
    
    ####################################################################
    # 2023-11: exit before making Incremental plots
    ####################################################################
    
    if True:
        return tracker.get_files()

    # Note: The code below is unreachable (as in the original MATLAB)
    # but is kept for reference in case incremental plots are needed later
    
    ####################################################################
    # Slews and Characterization vs. time -- Incremental
    ####################################################################
    
    fig, ax = plt.subplots(figsize=(8.5, 5))
    
    names = ['h_time_det_incr', 'h_time_char_incr', 'h_time_slew_incr']
    names_legend = ['Detection Time', 'Characterization Time', 'Slew Time']
    n_plot = len(names)
    
    # Put each of the above detection-times on one plot
    for n, name in enumerate(names):
        f_mean = f'{name}_mean'
        f_std = f'{name}_std'
        
        # Note: thinner linewidth for incremental plots
        ebar_props_incr = ebar_props.copy()
        ebar_props_incr['linewidth'] = 1
        
        ax.errorbar(tsamp + t_offsets[n],
                    t_det_time[f_mean].values,
                    yerr=t_det_time[f_std].values,
                    color=line_colors[n],
                    **ebar_props_incr)
    
    style_det_plot(ax,
                  'Mission Observation Scheduling (Incremental) vs. Mission Time',
                  'Incremental Time Used [days/month]',
                  names_legend)
    write_plots(fig, 'obstime-incr')
    plt.close(fig)


def main():
    """
    Command-line interface for plot_drm_time_used
    """
    parser = argparse.ArgumentParser(
        description='Plot time-used in a DRM-set',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
    python plot_drm_time_used.py "data/%s.%s" "output/det-%s.%s"
    
The first argument is the source template with two %%s placeholders that will be
filled with ("info", "csv") and ("det_time", "csv") to locate input files.

The second argument is the destination template with two %%s placeholders for
the plot name and file extension.
        """
    )
    
    parser.add_argument('src_tmpl', type=str,
                       help='Source template string (e.g., "data/%%s.%%s")')
    parser.add_argument('dest_tmpl', type=str,
                       help='Destination template string (e.g., "output/det-%%s.%%s")')
    parser.add_argument('--mode_op', type=str, default='*',
                       help='Operation mode (default: "*")')
    
    args = parser.parse_args()
    
    # Create mode dictionary
    mode = {'op': args.mode_op}
    
    # Run the plotting function
    try:
        plot_drm_time_used(args.src_tmpl, args.dest_tmpl, mode)
    except Exception as e:
        print(f"{PROGNAME}: Fatal: Unexpected error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
