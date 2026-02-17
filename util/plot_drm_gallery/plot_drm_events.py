#!/usr/bin/env python
"""
Plot event durations in a DRM-set

Plots of event-durations such as characterization integration times,
slews, etc.
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
except ImportError:
    import common_style as cs

# Program name for error messages
PROGNAME = os.path.basename(sys.argv[0])


def plot_drm_events(reduce_info, plot_data, dest_tmpl, mode):
    """
    Plot event durations in a drm-set

    Plots of event-durations such as characterization integration times,
    slews, etc.

    Parameters
    ----------
    reduce_info : dict
        Metadata dict from reduce-info.csv
    plot_data : list of DataFrame
        Pre-loaded CSV data [events]
    dest_tmpl : str
        Template string for output file paths with two %s placeholders
    mode : dict
        Dictionary with 'op' key containing operation mode string

    Outputs
    -------
    Saves plots to disk with names:
        duration-det-b1.png, duration-det-b2.png
        duration-char-b0.png, duration-char-b1.png, duration-char-b2.png
        duration-slew-b0.png, duration-slew-b1.png
    """

    # Unpack CSV data
    t_events, = plot_data
    
    # Track output files
    tracker = cs.PlotTracker(ext_list=mode.get('ext_list'))

    # Inner function: Set up plot/axis styles, title, axis labels
    def style_event_plot(ax, title1, xtext, ytext, legtext):
        """Style the event plot with title, labels, and legend"""
        # Clip Y range at 0
        ylim = ax.get_ylim()
        ax.set_ylim(max(0, ylim[0]), ylim[1])
        
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
        tracker.write_plots(fig, dest_name, dest_tmpl, verbose=mode['verbose'])

    # Sample times - three resolutions
    try:
        tsamp_0 = t_events['h_event_b0_duration_lo'].values
        tsamp_1 = t_events['h_event_b1_duration_lo'].values
        tsamp_2 = t_events['h_event_b2_duration_lo'].values
    except KeyError as e:
        print(f"{PROGNAME}: Must skip: Missing required column in events file: {e}", 
              file=sys.stderr)
        return []
    
    # Offsets on bars in the plot -- unused at present
    t_offsets_0 = [0]  # days units
    t_offsets_1 = [0]  # days units
    t_offsets_2 = [t / 24.0 for t in t_offsets_1]  # days - finer res. (?)
    
    # Note: the code below is set up for multiple-line error-bar
    # plots, but we have not used multiple bars.
    
    ####################################################################
    # Setup
    ####################################################################
    
    # Manual line color order
    # line_colors = ['tab:blue', 'tab:red', 'tab:orange']
    line_colors = ['tab:blue', 'tab:orange', 'black']

    # uniform error bar properties
    ebar_props = {"marker": '.',
                  "linewidth": 1.7,
                  "errorevery": 1, # highly variable: must do every one
                  "elinewidth": 1,
                  "ecolor": 'gray', # same ebar color ok for one-line plots
                  "capsize": 1.5}


    ####################################################################
    # Detection duration (daily binning)
    ####################################################################
    
    fig, ax = plt.subplots(figsize=(8.5, 5))
    
    # Just detection duration for now
    names = ['h_event_det_b1_duration']
    names_legend = ['Detection Time']
    n_plot = len(names)
    
    # Put each of the above detection-times on one plot
    for n, name in enumerate(names):
        f_mean = f'{name}_mean'
        f_std = f'{name}_std'
        
        ax.errorbar(tsamp_1 + t_offsets_1[n],
                   t_events[f_mean].values,
                   yerr=t_events[f_std].values,
                   color=line_colors[n],
                   **ebar_props)
    
    style_event_plot(ax, 'Mean Integration Time: Detection',
                    'Integration Time [day]',
                    'Frequency [density]', names_legend)
    write_plots(fig, 'duration-det-b1')
    plt.close(fig)
    
    ####################################################################
    # Detection duration (hourly binning)
    ####################################################################
    
    fig, ax = plt.subplots(figsize=(8.5, 5))
    
    # Scale factor (days -> hours)
    # We're plotting a density, so if we scale X, we have to scale Y
    # to preserve unit area under the density curve
    sf = 24
    
    names = ['h_event_det_b2_duration']
    names_legend = ['Detection Time']
    n_plot = len(names)
    
    # Put each of the above detection-times on one plot
    for n, name in enumerate(names):
        f_mean = f'{name}_mean'
        f_std = f'{name}_std'
        
        ax.errorbar(sf * (tsamp_2 + t_offsets_2[n]),
                   t_events[f_mean].values / sf,
                   yerr=t_events[f_std].values / sf,
                   color=line_colors[n],
                   **ebar_props)
    
    style_event_plot(ax, 'Mean Integration Time: Detection',
                    'Integration Time [hour]',
                    'Frequency [density]', names_legend)
    write_plots(fig, 'duration-det-b2')
    plt.close(fig)
    
    ####################################################################
    # Characterization duration (daily binning)
    ####################################################################
    
    fig, ax = plt.subplots(figsize=(8.5, 5))
    
    # Char duration
    names = ['h_event_char_b1_duration']
    names_legend = ['Characterization Time']
    n_plot = len(names)
    
    # Put each of the above time(s) on one plot
    for n, name in enumerate(names):
        f_mean = f'{name}_mean'
        f_std = f'{name}_std'
        
        ax.errorbar(tsamp_1 + t_offsets_1[n],
                   t_events[f_mean].values,
                   yerr=t_events[f_std].values,
                   color=line_colors[n],
                   **ebar_props)
    
    style_event_plot(ax, 'Mean Integration Time: Characterization',
                    'Integration Time [day]',
                    'Frequency [density]', names_legend)
    write_plots(fig, 'duration-char-b1')
    plt.close(fig)
    
    ####################################################################
    # Characterization duration (hourly binning)
    ####################################################################
    
    fig, ax = plt.subplots(figsize=(8.5, 5))
    
    # Scale factor (days -> hours)
    sf = 24
    
    names = ['h_event_char_b2_duration']
    names_legend = ['Characterization Time']
    n_plot = len(names)
    
    # Put each of the above times on one plot
    for n, name in enumerate(names):
        f_mean = f'{name}_mean'
        f_std = f'{name}_std'
        
        ax.errorbar(sf * (tsamp_2 + t_offsets_2[n]),
                   t_events[f_mean].values / sf,
                   yerr=t_events[f_std].values / sf,
                   color=line_colors[n],
                   **ebar_props)
    
    style_event_plot(ax, 'Mean Integration Time: Characterization',
                    'Integration Time [hour]',
                    'Frequency [density]', names_legend)
    write_plots(fig, 'duration-char-b2')
    plt.close(fig)
    
    ####################################################################
    # Characterization duration (2-day binning)
    ####################################################################
    
    fig, ax = plt.subplots(figsize=(8.5, 5))
    
    # Scale factor: plot is a density denominated in days,
    # but bin-width of x-axis is 2 days. This makes the
    # b1 plot and the b0 plot in the same units (density in days).
    sf = 2
    
    names = ['h_event_char_b0_duration']
    names_legend = ['Characterization Time']
    n_plot = len(names)
    
    # Put each of the above times on one plot
    for n, name in enumerate(names):
        f_mean = f'{name}_mean'
        f_std = f'{name}_std'
        
        ax.errorbar(1 * (tsamp_0 + t_offsets_0[n]),
                   sf * t_events[f_mean].values,
                   yerr=sf * t_events[f_std].values,
                   color=line_colors[n],
                   **ebar_props)
    
    style_event_plot(ax, 'Mean Integration Time: Characterization',
                    'Integration Time [day]',
                    'Frequency [density]', names_legend)
    write_plots(fig, 'duration-char-b0')
    plt.close(fig)
    
    ####################################################################
    # Slew duration (daily binning)
    ####################################################################
    
    fig, ax = plt.subplots(figsize=(8.5, 5))
    
    # Slew duration
    names = ['h_event_slew_b1_duration']
    names_legend = ['Slew Time']
    n_plot = len(names)
    
    # Put each of the above time(s) on one plot
    for n, name in enumerate(names):
        f_mean = f'{name}_mean'
        f_std = f'{name}_std'
        
        ax.errorbar(tsamp_1 + t_offsets_1[n],
                   t_events[f_mean].values,
                   yerr=t_events[f_std].values,
                   color=line_colors[n],
                   **ebar_props)
    
    style_event_plot(ax, 'Mean Times: Slew',
                    'Slew Time [day]',
                    'Frequency [density]', names_legend)
    write_plots(fig, 'duration-slew-b1')
    plt.close(fig)
    
    ####################################################################
    # Slew duration (2-day binning)
    ####################################################################
    
    fig, ax = plt.subplots(figsize=(8.5, 5))
    
    # Slew duration
    names = ['h_event_slew_b0_duration']
    names_legend = ['Slew Time']
    n_plot = len(names)
    
    # Scale factor: plot is a density denominated in days,
    # but bin-width of x-axis is 2 days. This makes the
    # b1 plot and the b0 plot in the same units (density in days).
    sf = 2
    
    # Put each of the above time(s) on one plot
    for n, name in enumerate(names):
        f_mean = f'{name}_mean'
        f_std = f'{name}_std'
        
        ax.errorbar(tsamp_0 + t_offsets_0[n],
                   sf * t_events[f_mean].values,
                   yerr=sf * t_events[f_std].values,
                   color=line_colors[n],
                   **ebar_props)
    
    style_event_plot(ax, 'Mean Times: Slew',
                    'Slew Time [day]',
                    'Frequency [density]', names_legend)
    write_plots(fig, 'duration-slew-b0')
    plt.close(fig)

    return tracker.get_files()


def main():
    """
    Command-line interface for plot_drm_events
    """
    parser = argparse.ArgumentParser(
        description='Plot event durations in a DRM-set',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
    python plot_drm_events.py "data/%s.%s"
    python plot_drm_events.py "data/%s.%s" "output/det-%s.%s"
    
The first argument is the source template with two %%s placeholders that will be
filled with ("info", "csv") and ("events", "csv") to locate input files.

The second argument is the destination template with two %%s placeholders for
the plot name and file extension (default: "scenario/gfx/det-%s.%s").
        """
    )
    
    parser.add_argument('src_tmpl', type=str,
                       help='Source template string (e.g., "data/%%s.%%s")')
    parser.add_argument('dest_tmpl', type=str, 
                       help='Destination template string (default: "scenario/gfx/det-%%s.%%s")')
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
    plot_data = cs.load_csv_files(args.src_tmpl, ['events'])
    rv = plot_drm_events(reduce_info, plot_data, args.dest_tmpl, mode)
    return rv


if __name__ == '__main__':
    rv = main()
    if rv is None:
        print(f"Plots failed. Error signaled.", file=sys.stderr)
    else:
        print(f"Done. Wrote {len(rv)} plot(s).")
    sys.exit(1 if rv is None else 0)
