#!/usr/bin/env python
"""
Plot yield, etc., against star targets in a drm-set

Make, and write out, plots of various quantities pertaining to each
star target in a set of DRMs.  

An example plot of this type is yield for each star, plotted against
stellar luminosity and distance.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import argparse
import sys
import os
# uniform plot appearance
from . import common_style as cs

# Program name for error messages
PROGNAME = os.path.basename(sys.argv[0])

# Verbosity (also set from mode)
VERBOSE = 1


def plot_drm_star_targets(src_tmpl, dest_tmpl, mode):
    """
    Plot yield, etc., against star targets in a drm-set
    
    Make, and write out, plots of various quantities pertaining to each
    star target in a set of DRMs.
    
    Parameters
    ----------
    src_tmpl : str
        Template string for input file paths with two %s placeholders
        Example: "data/%s.%s" will be filled with ("info", "csv") and ("star_targ", "csv")
    dest_tmpl : str
        Template string for output file paths with two %s placeholders
    mode : dict
        Dictionary with 'op' key containing operation mode string
        
    Outputs
    -------
    Saves plots to disk with names:
        perstar-det-allplan-cume.png, perstar-det-allplan-uniq.png
        perstar-det-earth-cume.png, perstar-det-earth-uniq.png
        perstar-char-allplan-cume.png, perstar-char-allplan-uniq.png
        perstar-char-earth-cume.png, perstar-char-earth-uniq.png
        perstar-det-t-int.png, perstar-char-t-int.png
        perstar-det-allplan-rank.png, perstar-det-allplan-frac.png
        perstar-char-allplan-rank.png, perstar-char-allplan-frac.png
        perstar-det-earth-rank.png, perstar-det-earth-frac.png
        perstar-char-earth-rank.png, perstar-char-earth-frac.png
        perstar-det-tint-yield.png, perstar-char-tint-yield.png
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
        star_targ_file = src_tmpl % ("star-target", "csv")
        t_star_targ = pd.read_csv(star_targ_file)
    except Exception as e:
        print('Star target plots: skipping (re-run reduction?).')
        return []

    # Allow skipping these plots so we don't fail on old reductions
    if t_star_targ.empty:
        print('Star target plots: skipping (re-run reduction?).')
        return []

    # Skip unless mode.op contains our name or a *
    if '*' not in mode.get('op', '') and 'perstar' not in mode.get('op', ''):
        print('Star target plots: skipping, as directed.')
        return []
    
    # Common data
    
    # Color map
    if False:
        # original, gives DeprecationWarning
        cmap_full = cm.get_cmap('jet')(np.linspace(0, 1, 280))
    else:
        cmap_full = plt.get_cmap('jet')(np.linspace(0, 1, 280))
    cmap = cmap_full[-256:, :]  # initial blue too dark
    cmap = cm.colors.ListedColormap(cmap)
    
    # 'scatter' dot size
    dot_size = 35
    
    # Clip target plots at 30 pc
    dist_axis_limit = [0.0, 30.0]
    
    # File extensions to write
    ext_list = ['png']

    # Track output files
    tracker = cs.PlotTracker()
    tracker.set_ext_list(ext_list)

    # Color for un-observed stars
    unseen_color = np.array([1, 1, 1]) * 0.7
    
    # Inner function: Set up plot/axis styles, title, axis labels
    def style_scatter_plot(ax, obs_name, obs_unit):
        """Style the scatter plot with title, labels, and colorbar"""
        # Plot/axis styles
        title1 = f'Luminosity vs. Distance, Shaded by {obs_name}'
        title2 = cs.plot_make_title(t_info)
        
        # Set title (preventing special interpretation of _) with bold
        ax.set_title(f'{title2}\n{title1}', fontsize=11*1.1, fontweight='bold')
        ax.set_xlabel('Distance [pc]', fontweight='bold')
        ax.set_ylabel('Luminosity [L$_{sun}$]', fontweight='bold')
        ax.tick_params(labelsize=13)
        ax.set_yscale('log')  # semilog
        ax.set_xlim(dist_axis_limit)  # clip distance
        ax.autoscale(enable=True, axis='y', tight=True)  # tight on Y
        # major grids
        ax.grid(True) 
        # minor y grid
        ax.grid(which='minor', axis='y', linestyle=':', linewidth='0.5', color='black', alpha=0.5)
        
        # Colorbar
        ylabel_string = f'{obs_name} [{obs_unit}]'
        cbar = plt.colorbar(ax.collections[-1], ax=ax)
        cbar.set_label(ylabel_string, fontweight='bold')
    
    # Inner function: write the current figure to files
    def write_plots(fig, file_tag):
        """Write the current figure to various files"""
        tracker.write_plots(fig, file_tag, dest_tmpl, verbose=VERBOSE)

    # Graphics setup
    
    # Zero-visit (unseen) stars are not the same as zero-yield stars
    # This variable keeps track of what stars were visited
    try:
        seen = (t_star_targ['h_star_det_visit_mean'].values > 0)
        earth = (t_star_targ['h_star_det_earth_cume_mean'].values > 0)
    except KeyError as e:
        print(f"{PROGNAME}: Fatal: Missing required column in star_targ file: {e}", 
              file=sys.stderr)
        sys.exit(1)
    
    # The plot_menu variable controls the properties that are plotted
    # in standard distance/luminosity coordinates. Each entry is one
    # plot, with a variable to plot, units, filename output, etc.
    #
    # Variable format in each row below:
    #   'h_star_...'    -> entry in the t_star_targ table
    #   passthru/np.log10 -> variable transform for the color-coding
    #   True/False -> put on the "we visited this" overlay (for Earths)
    #   'Title'    -> plot title
    #   'unit'     -> plot color-coding label (x/y axis are always dist/lum)
    #   'perstar_...' -> output filename component
    
    # Dummy function to pass matrix through unchanged
    passthru = lambda x: x
    
    plot_menu = [
        ['h_star_det_plan_cume_mean', passthru, False,
         'Mean Total Detections', 'count',
         'perstar-det-allplan-cume'],
        ['h_star_det_plan_uniq_mean', passthru, False,
         'Mean Unique Detections', 'count',
         'perstar-det-allplan-uniq'],
        ['h_star_det_earth_cume_mean', passthru, True,
         'Mean Total Detections: Earth', 'count',
         'perstar-det-earth-cume'],
        ['h_star_det_earth_uniq_mean', passthru, True,
         'Mean Unique Detections: Earth', 'count',
         'perstar-det-earth-uniq'],
        ['h_star_char_plan_cume_mean', passthru, False,
         'Mean Total Characterizations', 'count',
         'perstar-char-allplan-cume'],
        ['h_star_char_plan_uniq_mean', passthru, False,
         'Mean Unique Characterizations', 'count',
         'perstar-char-allplan-uniq'],
        ['h_star_char_earth_cume_mean', passthru, True,
         'Mean Total Characterizations: Earth', 'count',
         'perstar-char-earth-cume'],
        ['h_star_char_earth_uniq_mean', passthru, True,
         'Mean Unique Characterizations: Earth', 'count',
         'perstar-char-earth-uniq'],
        ['h_star_det_tInt_mean', passthru, False,
         'Mean Integration Time (Det.)', 'day',
         'perstar-det-t-int'],
        ['h_star_char_tInt_mean', passthru, False,
         'Mean Integration Time (Char.)', 'day',
         'perstar-char-t-int'],
        ['h_star_det_plan_value_mean', np.log10, False,
         'Detection Rank', 'log$_{10}$ count/day',
         'perstar-det-allplan-rank'],
        ['h_star_det_plan_frac_mean', passthru, False,
         'Planets Detected/Planets Present', 'count/count',
         'perstar-det-allplan-frac'],
        ['h_star_char_plan_value_mean', np.log10, False,
         'Characterization Rank', 'log$_{10}$ count/day',
         'perstar-char-allplan-rank'],
        ['h_star_char_plan_frac_mean', passthru, False,
         'Planets Characterized/Planets Present', 'count/count',
         'perstar-char-allplan-frac'],
        ['h_star_det_earth_value_mean', np.log10, True,
         'Earth Detection Rank', 'log$_{10}$ count/day',
         'perstar-det-earth-rank'],
        ['h_star_det_earth_frac_mean', passthru, True,
         'Earths Detected/Earths Present', 'count/count',
         'perstar-det-earth-frac'],
        ['h_star_char_earth_value_mean', np.log10, True,
         'Earth Characterization Rank', 'log$_{10}$ count/day',
         'perstar-char-earth-rank'],
        ['h_star_char_earth_frac_mean', passthru, True,
         'Earths Characterized/Earths Present', 'count/count',
         'perstar-char-earth-frac'],
    ]
    
    # This loop runs through the table above
    for plot_num, plot_spec in enumerate(plot_menu):
        if len(plot_spec) != 6:
            print(f"{PROGNAME}: Fatal: Plot selector {plot_num} has wrong number of elements",
                  file=sys.stderr)
            sys.exit(1)
        
        plot_prop = plot_spec[0]   # property to plot
        plot_xform = plot_spec[1]  # transformation of the property
        plot_earth = plot_spec[2]  # is it a plot of Earthlike dets/chars?
        plot_title = plot_spec[3]  # title, text
        plot_unit = plot_spec[4]   # unit of plot_prop above
        plot_file = plot_spec[5]   # filename
        
        fig, ax = plt.subplots(figsize=(8.5, 5))
        
        # Selector for the final, "on-top" scatter plot
        if plot_earth:
            selector = earth
        else:
            selector = seen
        
        # Un-visited in plain color, transparency, smaller dots
        ax.scatter(t_star_targ.loc[~seen, 'star_dist'].values,
                  t_star_targ.loc[~seen, 'star_L'].values,
                  s=dot_size/4,
                  c=[unseen_color],
                  marker='.',
                  alpha=0.2)
        
        # Show visited stars, Earth or otherwise
        if plot_earth:
            ax.scatter(t_star_targ.loc[seen, 'star_dist'].values,
                      t_star_targ.loc[seen, 'star_L'].values,
                      s=dot_size/2,
                      c=[[0.4, 0.4, 0.4]],
                      marker='.',
                      alpha=0.2)
        
        # Visited, in color varying by detection rate
        scatter = ax.scatter(t_star_targ.loc[selector, 'star_dist'].values,
                           t_star_targ.loc[selector, 'star_L'].values,
                           s=dot_size,
                           c=plot_xform(t_star_targ.loc[selector, plot_prop].values),
                           cmap=cmap,
                           marker='.')
        
        # Plot/axis styles
        style_scatter_plot(ax, plot_title, plot_unit)
        write_plots(fig, plot_file)
        plt.close(fig)
    
    ####################################################################
    # Change to another plot format
    ####################################################################
    
    # Alternative Plot: tInt vs. yield -- Detections
    ####################################################################
    
    fig, ax = plt.subplots(figsize=(8.5, 5))
    
    scatter = ax.scatter(t_star_targ.loc[seen, 'h_star_det_tInt_mean'].values,
                        t_star_targ.loc[seen, 'h_star_det_plan_cume_mean'].values,
                        s=dot_size,
                        c=t_star_targ.loc[seen, 'h_star_det_tobs1_mean'].values,
                        cmap=cmap,
                        marker='.')
    
    # log for tInt
    ax.set_xscale('log')
    
    # Plot/axis styles
    title1 = 'Cumulative Integration Time vs. Detections\nShaded by First-Observation Time'
    title2 = cs.plot_make_title(t_info)
    
    ax.set_title(f'{title2}\n{title1}', fontsize=11*1.1, fontweight='bold')
    ax.set_xlabel('Integration Time [d]', fontweight='bold')
    ax.set_ylabel('Mean Yield [count]', fontweight='bold')
    ax.tick_params(labelsize=13)
    ax.grid(True)
    # minor x grid (log [d])
    ax.grid(which='minor', axis='x', linestyle=':', linewidth='0.5', color='black', alpha=0.5)
    
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('Time of First Obs. [day]', fontweight='bold')
    
    write_plots(fig, 'perstar-det-tint-yield')
    plt.close(fig)
    
    # Alternative Plot: tInt vs. yield -- Characterizations
    ####################################################################
    
    fig, ax = plt.subplots(figsize=(8.5, 5))
    
    scatter = ax.scatter(t_star_targ.loc[seen, 'h_star_char_tInt_mean'].values,
                        t_star_targ.loc[seen, 'h_star_char_plan_cume_mean'].values,
                        s=dot_size,
                        c=t_star_targ.loc[seen, 'h_star_char_tobs1_mean'].values,
                        cmap=cmap,
                        marker='.')
    
    # log for tInt
    ax.set_xscale('log')
    
    # Plot/axis styles
    title1 = 'Cumulative Integration Time vs. Characterizations\nShaded by First-Observation Time'
    title2 = cs.plot_make_title(t_info)
    
    ax.set_title(f'{title2}\n{title1}', fontsize=11*1.1, fontweight='bold')
    ax.set_xlabel('Integration Time [d]', fontweight='bold')
    ax.set_ylabel('Mean Characterizations [count]', fontweight='bold')
    ax.tick_params(labelsize=13)
    ax.grid(True)
    # minor x grid (log [d])
    ax.grid(which='minor', axis='x', linestyle=':', linewidth='0.5', color='black', alpha=0.5)
    
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('Time of First Obs. [day]', fontweight='bold')
    
    write_plots(fig, 'perstar-char-tint-yield')
    plt.close(fig)

    return tracker.get_files()


def main():
    """
    Command-line interface for plot_drm_star_targets
    """
    parser = argparse.ArgumentParser(
        description='Plot yield, etc., against star targets in a DRM-set',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
    python plot_drm_star_targets.py "data/%s.%s" "output/det-%s.%s"
    
The first argument is the source template with two %%s placeholders that will be
filled with ("info", "csv") and ("star_targ", "csv") to locate input files.

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
        plot_drm_star_targets(args.src_tmpl, args.dest_tmpl, mode)
    except Exception as e:
        print(f"{PROGNAME}: Fatal: Unexpected error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
