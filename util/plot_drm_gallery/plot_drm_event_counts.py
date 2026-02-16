#!/usr/bin/env python
"""
Plot event counts in a drm-set

Plots of event-counts such as number of slews, characterizations, 
detections, averaged over an ensemble.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
import os
from scipy.ndimage import uniform_filter1d
# uniform plot appearance
from . import common_style as cs

# Program name for error messages
PROGNAME = os.path.basename(sys.argv[0])

# Verbosity (also set from mode)
VERBOSE = 1

def plot_drm_event_counts(reduce_info, src_tmpl, dest_tmpl, mode):
    """
    Plot event counts in a drm-set
    
    Plots of event-counts such as number of slews, characterizations, 
    detections, averaged over an ensemble.
    
    Parameters
    ----------
    src_tmpl : str
        Template string for input file paths with two %s placeholders
        Example: "data/%s.%s" will be filled with ("info", "csv"), 
        ("counts", "csv"), and ("earth_counts", "csv")
    dest_tmpl : str
        Template string for output file paths with two %s placeholders
    mode : dict
        Dictionary with 'op' key containing operation mode string
        
    Outputs
    -------
    Saves plots to disk with names:
        event-count-slew.png
        event-count-det.png, event-count-det-zoom.png
        event-count-det-rv.png
        earth-char-count.png
        earth-char-count-strict.png, earth-char-count-all.png
    """
    
    # Load data using the source template
    # update global verbosity
    global VERBOSE
    VERBOSE = mode.get('verbose', VERBOSE)

    try:
        counts_file = src_tmpl % ("event-counts", "csv")
        t_counts = pd.read_csv(counts_file)
    except Exception as e:
        print(f"{PROGNAME}: Fatal: Could not load counts file: {e}", file=sys.stderr)
        sys.exit(1)
    
    try:
        earth_counts_file = src_tmpl % ("earth-char-count", "csv")
        t_earth_counts = pd.read_csv(earth_counts_file)
    except Exception as e:
        print(f"Warning: Could not load earth-counts file: {e}", file=sys.stderr)
        t_earth_counts = pd.DataFrame()
    
    # Skip unless mode.op contains our name or a *
    if '*' not in mode.get('op', '') and 'event-counts' not in mode.get('op', ''):
        print('Event-count plots: skipping, as directed.')
        return []
    
    # File extensions to write
    ext_list = ['png']

    # Track output files
    tracker = cs.PlotTracker()
    tracker.set_ext_list(ext_list)

    # Inner function: Set up plot/axis styles, title, axis labels
    def style_count_plot(ax, title1, xtext, ytext, legtext):
        """Style the count plot with title, labels, and legend"""
        # Clip Y range at 0
        ylim = ax.get_ylim()
        ax.set_ylim(max(0, ylim[0]), ylim[1])
        
        # Clip X at 0 - was seeing significant negative extension
        xlim = ax.get_xlim()
        ax.set_xlim(0, xlim[1])
        # Format the title
        title2 = cs.plot_make_title(reduce_info)
        
        # Set title (preventing special interpretation of _) with bold
        ax.set_title(f'{title2}\n{title1}', fontsize=11*1.1, fontweight='bold')
        
        # Axis labels
        ax.set_xlabel(xtext, fontweight='bold')
        ax.set_ylabel(ytext, fontweight='bold')
        
        # NE corner is good for histograms
        if len(legtext) > 0:
            ax.legend(legtext, loc='upper right')
        
        ax.tick_params(labelsize=13)
        ax.grid(True)
    
    # Inner function: write the current figure to files
    def write_plots(fig, dest_name):
        """Write the current figure to various files"""
        tracker.write_plots(fig, dest_name, dest_tmpl, verbose=VERBOSE)

    # Histogram domain (in counts)
    try:
        ct_samp_1 = t_counts['h_event_count_lo'].values
    except KeyError as e:
        print(f"{PROGNAME}: Fatal: Missing required column in counts file: {e}", 
              file=sys.stderr)
        sys.exit(1)
    
    # Offsets on bars in the plot -- unused at present
    ct_offsets_1 = [0]  # counts units
    
    # Manual line color order
    # line_colors = ['tab:blue', 'tab:red', 'tab:orange']
    line_colors = ['tab:blue', 'tab:green', 'tab:orange']

    # uniform error bar properties
    plot_props = {"marker": '.',
                  "linewidth": 1.7}

    # Note: the code below is set up for multiple-line error-bar
    # plots, but we have not used multiple bars.
    
    ####################################################################
    # Slews
    ####################################################################
    
    fig, ax = plt.subplots(figsize=(8.5, 5))
    
    # Just slew count in this plot
    names = ['h_event_count_slew']
    names_legend = ['Slews']
    n_plot = len(names)
    # Indexes to plot
    inx = slice(0, 150)
    
    # only using one for now
    slew_colors = ['tab:purple', 'tab:gray']

    # Put the above-selected counts on one plot
    for n, name in enumerate(names):
        f_mean = f'{name}_mean'
        f_std = f'{name}_std'
        
        extra_props = {'zorder': 2} if 'char' in name else {'zorder': 1}

        ax.plot(ct_samp_1[inx] + ct_offsets_1[n],
                t_counts[f_mean].values[inx],
                color=slew_colors[n],
                **plot_props)
        # (error bars just confuse this plot)
    
    style_count_plot(ax, 'Mean Event Count: Slews',
                    'Number of Events [count]',
                    'Frequency [density]', names_legend)
    write_plots(fig, 'event-count-slew')
    plt.close(fig)
    
    ####################################################################
    # Detections
    ####################################################################
    
    fig, ax = plt.subplots(figsize=(8.5, 5))
    
    # Detections, etc.
    names = ['h_event_count_det', 'h_event_count_char', 'h_event_count_detp']
    names_legend = ['Detections', 'Characterizations', 'Detections, Promoted']
    n_plot = len(names)
    # Indexes to plot (event counts)
    inx = slice(0, 500)
    
    # Put the above-selected counts on one plot
    for n, name in enumerate(names):
        f_mean = f'{name}_mean'
        f_std = f'{name}_std'
        
        # Average in bins, standard deviations average in quadrature
        # Using uniform_filter1d for moving average (equivalent to movmean)
        mean_avg = uniform_filter1d(t_counts[f_mean].values, size=3, mode='nearest')
        std_avg = np.sqrt(np.maximum(0.0, 
            uniform_filter1d(t_counts[f_std].values ** 2, size=3, mode='nearest')))
        
        extra_props = {'zorder': 2} if 'char' in name else {'zorder': 1}
        ax.plot(ct_samp_1[inx], mean_avg[inx],
                color=line_colors[n],
                **extra_props, 
                **plot_props)
        # (skip the error bars)
    
    style_count_plot(ax, 'Mean Event Count: Detection and Characterization Attempts',
                    'Number of Events [count]',
                    'Frequency [density]', names_legend)
    write_plots(fig, 'event-count-det')
    
    # Zoomed version
    style_count_plot(ax, 'Mean Event Count: Detection and Characterization Attempts (Zoomed)',
                    'Number of Events [count]',
                    'Frequency [density]', names_legend)
    ax.set_xlim(left=None, right=100)
    ax.set_ylim(bottom=0, top=None)
    write_plots(fig, 'event-count-det-zoom')
    plt.close(fig)
    
    ####################################################################
    # Detections - RV precursors
    ####################################################################
    
    fig, ax = plt.subplots(figsize=(8.5, 5))
    
    # Dets/chars -- plot order changed in python vs. matlab to harmonize "char" colors
    names = ['h_event_count_det_rvplan', 'h_event_count_char_rvplan', 'h_event_count_char']
    names_legend = ['Detections, RV Targets', 'Characterizations, RV Targets', 'Characterizations, Any']
    n_plot = len(names)
    # Indexes to plot (event counts)
    inx = slice(0, 500)
    
    # Ordinary chars are plotted last, and take color index 1
    # special color for RV chars
    # FIXME: not working RN
    rv_line_colors = [line_colors[0],  'tab:olive', line_colors[1]]

    # This guard is needed for old reductions
    if 'h_event_count_det_rvplan_mean' in t_counts.columns:
        
        # Put the above-selected counts on one plot
        for n, name in enumerate(names):
            f_mean = f'{name}_mean'
            f_std = f'{name}_std'
            
            # Average in bins, standard deviations average in quadrature
            mean_avg = uniform_filter1d(t_counts[f_mean].values, size=3, mode='nearest')
            std_avg = np.sqrt(np.maximum(0.0, 
                uniform_filter1d(t_counts[f_std].values ** 2, size=3, mode='nearest')))
            
            # chars on top, 
            extra_props = {
                'zorder': 2 if 'char' in name else 1,
                'linestyle': ':' if 'char_rv' in name else '-'
                }
            ax.plot(ct_samp_1[inx], mean_avg[inx],
                    color=line_colors[n],
                    **extra_props,
                    **plot_props)
            # (skip the error bars)
        
        style_count_plot(ax, 'Mean Event Count: RV Detection and Characterization Attempts',
                        'Number of Events [count]',
                        'Frequency [density]', names_legend)
        ax.set_xlim(left=None, right=30)
        ax.set_ylim(bottom=0, top=None)
        write_plots(fig, 'event-count-det-rv')
    
    plt.close(fig)
    
    ####################################################################
    # Earth Characterization Counts
    ####################################################################
    
    # Allow early exit if data not present
    if t_earth_counts.empty:
        return tracker.get_files()
    
    fig, ax = plt.subplots(figsize=(8.5, 5))
    
    # Histogram domain (in counts)
    # Bins include lower boundary, but exclude upper boundary
    # Thus: if "lo" is [0,1,2,...] then first bin is [0,1).
    try:
        ct_samp_1 = t_earth_counts['h_earth_char_count_lo'].values
    except KeyError as e:
        print(f"{PROGNAME}: Fatal: Missing required column in earth_counts file: {e}", 
              file=sys.stderr)
        sys.exit(1)
    
    # Offsets on bars in the plot -- unused at present
    ct_offsets_1 = [0, 0]  # counts units
    
    # Earth char counts appearing in this plot
    names = ['h_earth_char_all', 'h_earth_char_strict']
    names_legend = ['Earths Characterized (full/partial, any band)',
                    'Earths Characterized (full, all bands)']
    n_plot = len(names)
    
    # Put the above-selected counts on one plot
    for n, name in enumerate(names):
        f_mean = f'{name}_mean'
        f_std = f'{name}_std'
        
        # Guard against out-of-date csv files for strict mode
        if f_mean not in t_earth_counts.columns:
            print(f'Skipping earth chars ({f_mean}): redo make reduce to fix.')
            continue
        
        ax.plot(ct_samp_1 + ct_offsets_1[n],
               t_earth_counts[f_mean].values,
               color=line_colors[n],
                **plot_props)
        # (error bars just confuse this plot)
    
    # Clip the x-range, ensure 0 gets an explicit bin
    ax.set_xlim(-0.5, 50)
    ax.set_ylim(bottom=0, top=None)
    style_count_plot(ax, 'Number of Earths Characterized',
                    'Number of Earths [count]',
                    'Frequency [density]', names_legend)
    write_plots(fig, 'earth-char-count')
    plt.close(fig)
    
    ####################################################################
    # Version 2/a: Showing one mode only, bars (strict)
    ####################################################################
    
    fig, ax = plt.subplots(figsize=(8.5, 5))
    
    bar_color = [70/255, 130/255, 180/255]  # SteelBlue
    
    # Earth char counts appearing in this plot
    names = ['h_earth_char_strict']
    names_legend = ''  # 1 set of bars -> suppress legend
    n_plot = len(names)
    
    # Put the above-selected counts on one plot
    # *** NB: this bar-plot will fail for N_plot > 1
    for n, name in enumerate(names):
        f_mean = f'{name}_mean'
        f_std = f'{name}_std'
        
        # Guard against out-of-date csv files for strict mode
        if f_mean not in t_earth_counts.columns:
            print(f'Skipping earth chars ({f_mean}): redo make reduce to fix')
            continue
        
        ax.bar(ct_samp_1 + ct_offsets_1[n],
              t_earth_counts[f_mean].values,
              color=bar_color)
        # (error bars just confuse this plot)
    
    # Clip the x-range
    ax.set_xlim(-0.5, 50)
    ax.set_ylim(bottom=0, top=None)
    style_count_plot(ax, 'Number of Earths Characterized',
                    'Number of Earths [count]',
                    'Frequency [density]', names_legend)
    write_plots(fig, 'earth-char-count-strict')
    plt.close(fig)
    
    ####################################################################
    # Version 2/b: Showing one mode only, bars (all)
    ####################################################################
    
    fig, ax = plt.subplots(figsize=(8.5, 5))
    
    # Earth char counts appearing in this plot
    names = ['h_earth_char_all']
    names_legend = ''  # 1 set of bars -> suppress legend
    n_plot = len(names)
    
    # Put the above-selected counts on one plot
    # *** NB: this bar-plot will fail for N_plot > 1
    for n, name in enumerate(names):
        f_mean = f'{name}_mean'
        f_std = f'{name}_std'
        
        # Guard against out-of-date csv files for strict mode
        if f_mean not in t_earth_counts.columns:
            print(f'Skipping earth chars ({f_mean}): redo make reduce to fix')
            continue
        
        ax.bar(ct_samp_1 + ct_offsets_1[n],
              t_earth_counts[f_mean].values,
              color=bar_color)
        # (error bars just confuse this plot)
    
    # Clip the x-range
    ax.set_xlim(-0.5, 50)
    ax.set_ylim(bottom=0, top=None)
    style_count_plot(ax, 'Number of Earths Characterized',
                    'Number of Earths [count]',
                    'Frequency [density]', names_legend)
    write_plots(fig, 'earth-char-count-all')
    plt.close(fig)

    return tracker.get_files()


def main():
    """
    Command-line interface for plot_drm_event_counts
    """
    parser = argparse.ArgumentParser(
        description='Plot event counts in a DRM-set',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
    python plot_drm_event_counts.py "data/%s.%s" "output/det-%s.%s"
    
The first argument is the source template with two %%s placeholders that will be
filled with ("info", "csv"), ("counts", "csv"), and ("earth_counts", "csv") to 
locate input files.

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

    # Read info file and convert to dict
    info_file = args.src_tmpl % ("info", "csv")
    reduce_info = pd.read_csv(info_file).iloc[0].to_dict()

    # Run the plotting function
    try:
        plot_drm_event_counts(reduce_info, args.src_tmpl, args.dest_tmpl, mode)
    except Exception as e:
        print(f"{PROGNAME}: Fatal: Unexpected error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
