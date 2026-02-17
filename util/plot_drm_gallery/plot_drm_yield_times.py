#!/usr/bin/env python
"""
Plot detection/characterization times in a drm-set

Time-series plots of detections, characterizations - plotted vs. 
mission elapsed time.

Note: This plot family originally included only detections (cume,
uniq, revisit) vs. time - within the t_det_time table. We later created 
a newer table (t_yield_time) that has det/char, and allplanets/earth,
and multiband char info. This latter table is the one we're using.
That is, for these plots, we use reduce-yield-time.csv 
but *not* reduce-times.csv.

Produces plots entitled: 
 det-time-char-allplan-full-cume.png
 det-time-char-allplan-part-cume.png
 det-time-det-allplan-cume.png
 det-time-det-earth-cume.png
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


def plot_drm_yield_times(reduce_info, plot_data, dest_tmpl, mode):
    """
    Plot detection/characterization times in a drm-set

    Time-series plots of detections, characterizations - plotted vs.
    mission elapsed time.

    Parameters
    ----------
    reduce_info : dict
        Metadata dict from reduce-info.csv
    plot_data : list of DataFrame
        Pre-loaded CSV data [yield-time]
    dest_tmpl : str
        Template string for output file paths with two %s placeholders
    mode : dict
        Dictionary with 'op' key containing operation mode string

    Outputs
    -------
    Saves plots to disk with names:
        det-time-char-allplan-full-cume.png
        det-time-char-allplan-part-cume.png
        det-time-det-allplan-cume.png
        det-time-det-earth-cume.png
    """

    # Make extra plots?
    do_incremental_plot = '+' in mode.get('op', '')
    
    # Unpack CSV data
    t_yield_time, = plot_data
    
    # File extensions to write
    ext_list = ['png']

    # Track output files
    tracker = cs.PlotTracker()
    tracker.set_ext_list(ext_list)

    # Inner function: Set up plot/axis styles, title, axis labels
    def style_yield_plot(ax, title1, ytext, legtext):
        """Style the yield plot with title, labels, and legend"""
        # Clip Y range at 0
        ylim = ax.get_ylim()
        ax.set_ylim(max(0, ylim[0]), ylim[1])
        
        # tight x axis limits - was losing too much space
        ax.autoscale(enable=True, axis='x', tight=True)

        # Format the title
        title2 = cs.plot_make_title(reduce_info)
        
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
        tracker.write_plots(fig, dest_name, dest_tmpl, verbose=mode['verbose'])

    # Offsets on various error-bars in the plot 
    # (units of days with one complete sample per ~30 days)
    t_offsets = [0, 10, 5]
    
    # Sample times
    try:
        tsamp = t_yield_time['h_det_time_lo'].values
    except KeyError as e:
        print(f"{PROGNAME}: Will skip: Missing required column in yield_time file: {e}", 
              file=sys.stderr)
        return []
    
    # Manual line color order
    # line_colors = ['tab:blue', 'tab:red', 'tab:orange']
    line_colors = ['tab:blue', 'tab:green', 'tab:orange']

    # Uniform error bar properties
    ebar_props = {"marker": '.',
                  "linewidth": 1.7,
                  "errorevery": 1,
                  "elinewidth": 1,
                  "capsize": 1.5}
    
    ####################################################################
    # Detections vs. time
    ####################################################################
    
    # Needed for each item here: 
    # {fieldnames, plot title, filename}
    
    # NOTE: 
    # * The third plot entry is *only full chars*.  
    # * The last two plot entries *include full and partial chars*.  
    # I.e., both full and partial chars are included in the latter two.  
    # This extra detail in the plot title was removed to
    # make the plots in the SDET report less overwhelming.
    
    fig_roster = [
        {
            'names': ['h_time_det_allplan_cume', 'h_time_det_allplan_uniq', 
                     'h_time_det_allplan_revi'],
            'title': 'Detections [All]',
            'fname': 'time-det-allplan'
        },
        {
            'names': ['h_time_det_earth_cume', 'h_time_det_earth_uniq', 
                     'h_time_det_earth_revi'],
            'title': 'Detections [Earths]',
            'fname': 'time-det-earth'
        },
        {
            'names': ['h_time_char_full_allplan_cume_union',
                     'h_time_char_full_allplan_uniq_union', 
                     'h_time_char_full_allplan_revi_union'],
            'title': 'Characterizations [All Planets]',
            'fname': 'time-char-allplan-full'
        },
        {
            'names': ['h_time_char_part_allplan_cume_union',
                     'h_time_char_part_allplan_uniq_union', 
                     'h_time_char_part_allplan_revi_union'],
            'title': 'Characterizations [All Planets]',
            'fname': 'time-char-allplan-part'
        },
        {
            'names': ['h_time_char_part_earth_cume_union',
                     'h_time_char_part_earth_uniq_union', 
                     'h_time_char_part_earth_revi_union'],
            'title': 'Characterizations [Earths]',
            'fname': 'time-char-earth-part'
        }
    ]
    
    # 2023-11: typically skip the monthly a/k/a incremental plots
    
    # Iterate over all figures in the roster
    for fig_info in fig_roster:
        names = fig_info['names']
        dtxt = fig_info['title']
        fname = fig_info['fname']
        
        n_plot = len(names)
        
        # Skip unless the required fieldnames are present
        # (e.g., no char summaries are in t_yield_time for det-only missions)
        skipping = False
        for name in names:
            f_mean = f'{name}_mean'
            if f_mean not in t_yield_time.columns:
                skipping = True
                break
        
        if skipping:
            print(f'\tSkipping {fname} plot (missing field)')
            continue
        
        # Legend names
        names_legend = [
            f'All {dtxt}',
            f'Unique {dtxt}',
            f'Revisit {dtxt}'
        ]
        
        # A: Monthly plot (currently disabled)
        if do_incremental_plot:
            fig, ax = plt.subplots(figsize=(8.5, 5))
            
            # Plot each timeseries
            for n, name in enumerate(names):
                f_mean = f'{name}_mean'
                f_std = f'{name}_std'
                
                ax.errorbar(tsamp + t_offsets[n], 
                            t_yield_time[f_mean].values,
                            yerr=t_yield_time[f_std].values,
                            color=line_colors[n],
                            **ebar_props)
            
            style_yield_plot(
                ax,
                f'Monthly {dtxt} vs. Mission Time',
                f'{dtxt} [count/month]',
                names_legend
            )
            write_plots(fig, f'{fname}-month')
            plt.close(fig)
        
        # B: Cumulative plot
        # Note: we derive the error bars here
        
        fig, ax = plt.subplots(figsize=(8.5, 5))
        
        # Plot each timeseries with cumulative values
        for n, name in enumerate(names):
            f_mean = f'{name}_mean'
            f_std = f'{name}_std'
            
            # Cumulative mean, plus cumulative std error (adding in quadrature,
            # i.e., sqrt-sum-of-squares)
            cume_mean = np.cumsum(t_yield_time[f_mean].values)
            cume_std = np.sqrt(np.cumsum(t_yield_time[f_std].values ** 2))
            
            ax.errorbar(tsamp + t_offsets[n],
                        cume_mean,
                        yerr=cume_std,
                        color=line_colors[n],
                        **ebar_props)
        
        style_yield_plot(
            ax,
            f'Cumulative {dtxt} vs. Mission Time',
            f'{dtxt} [count]',
            names_legend
        )
        write_plots(fig, f'{fname}-cume')
        plt.close(fig)

    return tracker.get_files()


def main():
    """
    Command-line interface for plot_drm_yield_times
    """
    parser = argparse.ArgumentParser(
        description='Plot detection/characterization times in a DRM-set',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
    python plot_drm_yield_times.py "data/%s.%s" "output/det-%s.%s"
    
The first argument is the source template with two %%s placeholders that will be
filled with ("info", "csv") and ("yield_time", "csv") to locate input files.

The second argument is the destination template with two %%s placeholders for
the plot name and file extension.
        """
    )
    
    parser.add_argument('src_tmpl', type=str,
                       help='Source template string (e.g., "data/%%s.%%s")')
    parser.add_argument('dest_tmpl', type=str,
                       help='Destination template string (e.g., "output/det-%%s.%%s")')
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
    plot_data = cs.load_csv_files(args.src_tmpl, ['yield-time'])
    rv = plot_drm_yield_times(reduce_info, plot_data, args.dest_tmpl, mode)
    return rv


if __name__ == '__main__':
    rv = main()
    if rv is None:
        print(f"Plots failed. Error signaled.", file=sys.stderr)
    else:
        print(f"Done. Wrote {len(rv)} plot(s).")
    sys.exit(1 if rv is None else 0)
