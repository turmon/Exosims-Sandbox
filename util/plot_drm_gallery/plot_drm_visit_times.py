#!/usr/bin/env python
"""
Plot detection/characterization visits in a drm-set

Time-series plots of visits for the purpose of detection or
characterization - (#visits) plotted vs. mission elapsed time.

Note: This plot family is similar to yield, but is *not* yield:
  t_yield_time -- yield (#planets) versus mission elapsed time. 
  t_visit_time -- star-visits versus mission elapsed time.
So, "visits" are counted whether they result in success or not, and 
visits do not count planets, they count periods-of-integration.
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


def plot_drm_visit_times(reduce_info, plot_data, dest_tmpl, mode):
    """
    Plot detection/characterization visits in a drm-set

    Time-series plots of visits for the purpose of detection or
    characterization - (#visits) plotted vs. mission elapsed time.

    Parameters
    ----------
    reduce_info : dict
        Metadata dict from reduce-info.csv
    plot_data : list of DataFrame
        Pre-loaded CSV data [visit-time]
    dest_tmpl : str
        Template string for output file paths with two %s placeholders
    mode : dict
        Dictionary with 'op' key containing operation mode string

    Outputs
    -------
    Saves plots to disk with names:
        visit-time-det-cume.png
        visit-time-char-cume.png
    """

    # Unpack CSV data
    t_visit_time, = plot_data
    
    # Allow skipping this way
    if '0' in mode.get('op', ''):
        print('Visit time plots: skipping, as directed.')
        return []
    
    # File extensions to write
    ext_list = ['png']

    # Track output files
    tracker = cs.PlotTracker()
    tracker.set_ext_list(ext_list)

    # Inner function: Set up plot/axis styles, title, axis labels
    def style_visit_plot(ax, title1, ytext, legtext):
        """Style the visit plot with title, labels, and legend"""
        # Clip Y range at 0
        ylim = ax.get_ylim()
        ax.set_ylim(max(0, ylim[0]), ylim[1])
        
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
        tsamp = t_visit_time['h_det_time_lo'].values
    except KeyError as e:
        print(f"{PROGNAME}: Fatal: Missing required column in visit_time file: {e}", 
              file=sys.stderr)
        sys.exit(1)
    
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
    # Plot roster
    ####################################################################
    
    # Needed for each item here: 
    # {fieldnames, plot title, filename}
    
    fig_roster = [
        {
            'names': ['h_visit_det_visit', 'h_visit_det_uniq', 'h_visit_det_revi'],
            'title': 'Detection-Mode Target Visits',
            'fname': 'visit-time-det'
        },
        {
            'names': ['h_visit_char_visit', 'h_visit_char_uniq', 'h_visit_char_revi'],
            'title': 'Characterization-Mode Target Visits',
            'fname': 'visit-time-char'
        }
    ]
    
    # Do not make incremental/monthly plots
    do_incremental_plot = False
    
    # Iterate over all figures in the roster
    for fig_info in fig_roster:
        names = fig_info['names']
        dtxt = fig_info['title']
        fname = fig_info['fname']
        
        n_plot = len(names)
        
        # Skip unless the required fieldnames are present
        # (e.g., no char summaries are in t_visit_time for det-only missions)
        skipping = False
        for name in names:
            f_mean = f'{name}_mean'
            if f_mean not in t_visit_time.columns:
                skipping = True
                break
        
        if skipping:
            print(f'\tSkipping {fname} plot (missing field)')
            continue
        
        # Legend names
        names_legend = [
            f'All {dtxt}',
            f'First-time {dtxt}',
            f'Return {dtxt}'
        ]
        
        # A: Monthly plot (currently disabled)
        if do_incremental_plot:
            fig, ax = plt.subplots(figsize=(8.5, 5))
            
            # Plot each timeseries
            for n, name in enumerate(names):
                f_mean = f'{name}_mean'
                f_std = f'{name}_std'
                
                ax.errorbar(tsamp + t_offsets[n], 
                            t_visit_time[f_mean].values,
                            yerr=t_visit_time[f_std].values,
                            color=line_colors[n],
                            **ebar_props)  # NB: skinny
            
            style_visit_plot(
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
            cume_mean = np.cumsum(t_visit_time[f_mean].values)
            cume_std = np.sqrt(np.cumsum(t_visit_time[f_std].values ** 2))
            
            ax.errorbar(tsamp + t_offsets[n],
                        cume_mean,
                        yerr=cume_std,
                        color=line_colors[n],
                        **ebar_props)
        
        style_visit_plot(
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
    Command-line interface for plot_drm_visit_times
    """
    parser = argparse.ArgumentParser(
        description='Plot detection/characterization visits in a DRM-set',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
    python plot_drm_visit_times.py "data/%s.%s" "output/det-%s.%s"
    
The first argument is the source template with two %%s placeholders that will be
filled with ("info", "csv") and ("visit_time", "csv") to locate input files.

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
    plot_data = cs.load_csv_files(args.src_tmpl, ['visit-time'])
    rv = plot_drm_visit_times(reduce_info, plot_data, args.dest_tmpl, mode)
    return rv


if __name__ == '__main__':
    rv = main()
    if rv is None:
        print(f"Plots failed. Error signaled.", file=sys.stderr)
    else:
        print(f"Done. Wrote {len(rv)} plot(s).")
    sys.exit(1 if rv is None else 0)
