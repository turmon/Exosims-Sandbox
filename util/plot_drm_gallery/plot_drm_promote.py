#!/usr/bin/env python
"""
Plot target promotion counts in a drm-set

Time-series plots of target promotion criteria
Histograms of planet-counts meeting criteria
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


def plot_drm_promote(reduce_info, plot_data, dest_tmpl, mode):
    """
    Plot target promotion counts in a drm-set

    Time-series plots of target promotion criteria
    Histograms of planet-counts meeting criteria

    Parameters
    ----------
    reduce_info : dict
        Metadata dict from reduce-info.csv
    plot_data : list of DataFrame
        Pre-loaded CSV data [promote, promote-hist]
    dest_tmpl : str
        Template string for output file paths with two %s placeholders
    mode : dict
        Dictionary with 'op' key containing operation mode string

    Outputs
    -------
    Saves plots to disk with names:
        promote-allplan-cume.png, promote-hzone-cume.png
        promote-earth-cume.png, promote-star-cume.png
        promote-star-span-cume.png
        phist-hzone-2year.png, phist-hzone-3year.png
        phist-earth-2year.png, phist-earth-3year.png
        phist-star-2year.png, phist-star-3year.png
        phist-star-span-2year.png, phist-star-span-3year.png
    """

    # Unpack CSV data
    t_promote, t_phist = plot_data

    # Skip plots if the input was empty
    if t_promote.empty or t_phist.empty:
        print('Promotion plots: No data. Skipping.')
        return []

    # Allow skipping this way
    if '0' in mode.get('op', ''):
        print('Promotion plots: skipping, as directed.')
        return []
    
    # File extensions to write
    ext_list = ['png']

    # Track output files
    tracker = cs.PlotTracker()
    tracker.set_ext_list(ext_list)

    # Inner function: Set up plot/axis styles, title, axis labels
    def style_promote_plot(ax, title1, ytext, legtext):
        """Style the promotion plot with title, labels, and legend"""
        # Clip Y range at 0
        ylim = ax.get_ylim()
        ax.set_ylim(max(0, ylim[0]), ylim[1])
        
        # Format the title
        title2 = cs.plot_make_title(reduce_info)
        
        # Set title (preventing special interpretation of _) with bold
        ax.set_title(f'{title2}\n{title1}', fontsize=11*1.1, fontweight='bold')
        
        # Axis labels
        ax.set_xlabel('Detection Integration Time [day]', fontweight='bold')
        ax.set_ylabel(ytext, fontweight='bold')
        
        # Legend
        ax.legend(legtext, loc='upper left')
        ax.tick_params(labelsize=13)
        ax.grid(True)
    
    # Inner function: Set up plot/axis styles for histogram plots
    def style_phist_plot(ax, title1, legtext):
        """Style the planet histogram plot with title, labels, and legend"""
        # Clip Y range at 0
        ylim = ax.get_ylim()
        ax.set_ylim(max(0, ylim[0]), ylim[1])
        
        # Format the title
        title2 = cs.plot_make_title(reduce_info)
        
        # Set title (preventing special interpretation of _) with bold
        ax.set_title(f'{title2}\n{title1}', fontsize=11*1.1, fontweight='bold')
        
        # Axis labels
        ax.set_xlabel('Planet Count [count]', fontweight='bold')
        ax.set_ylabel('Density', fontweight='bold')
        
        # Legend
        ax.legend(legtext, loc='upper left')
        ax.tick_params(labelsize=13)
        
        # Set marker and line properties
        for line in ax.get_lines():
            line.set_markersize(8)
            line.set_linewidth(2)
        
        ax.grid(True)
    
    # Inner function: write the current figure to files
    def write_plots(fig, dest_name):
        """Write the current figure to various files"""
        tracker.write_plots(fig, dest_name, dest_tmpl, verbose=mode['verbose'])

    # Offsets on various error-bars in the plot 
    # (units of days with one complete sample per ~30 days)
    # (fourth offset is rarely used)
    t_offsets = [0, 10, 5, -5]
    
    # Sample times
    try:
        tsamp = t_promote['h_promo_time_lo'].values
    except KeyError as e:
        print(f"{PROGNAME}: Fatal: Missing required column in promote file: {e}", 
              file=sys.stderr)
        sys.exit(1)
    
    # Manual line color order
    # line_colors = ['tab:blue', 'tab:red', 'tab:orange']
    line_colors = ['tab:blue', 'tab:green', 'tab:orange']
    
    # Uniform error bar properties
    ebar_props = {"marker": '.',
                  "linewidth": 1.7,
                  "errorevery": 4,
                  "elinewidth": 1,
                  "capsize": 1.5}
    
    ####################################################################
    # Promotions vs. time
    ####################################################################
    
    # A: Cumulative AllPlan plot
    ####################################################################
    
    fig, ax = plt.subplots(figsize=(8.5, 5))
    
    names = ['h_promo_count_allplan_cume', 'h_promo_span_allplan_cume', 'h_promo_promo_allplan_cume']
    names_legend = ['Obs. Count', 'Obs. Span', 'Promotion']
    n_plot = len(names)
    
    # Put each of the above detection-times on one plot
    skipping = False
    for n, name in enumerate(names):
        f_mean = f'{name}_q50'
        f_bar1 = f'{name}_q25'
        f_bar2 = f'{name}_q75'
        
        if f_mean not in t_promote.columns:
            # Stop at any error here
            print(f'{PROGNAME}: No {f_mean} in promotion table, skipping')
            skipping = True
            break
        
        # Asymmetric error bars
        yerr_lower = t_promote[f_mean].values - t_promote[f_bar1].values
        yerr_upper = t_promote[f_bar2].values - t_promote[f_mean].values
        
        ax.errorbar(tsamp + t_offsets[n],
                   t_promote[f_mean].values,
                   yerr=[yerr_lower, yerr_upper],
                   color=line_colors[n],
                   **ebar_props)  # NB: skinny
    
    # This happens when there are no DRM's, and in this case all fields are missing,
    # so all plots will fail
    if skipping:
        print(f'\t{PROGNAME}: Skipping promotion plots (re-run reduction, or no DRMs?).')
        return tracker.get_files()
    
    style_promote_plot(ax, 
                      'All Planets: Cumulative Promotions vs. Detector Time',
                      'Targets Passing Criterion [count]', 
                      names_legend)
    write_plots(fig, 'promote-allplan-cume')
    plt.close(fig)
    
    # B: Cumulative Hab Zone plot
    ####################################################################
    
    fig, ax = plt.subplots(figsize=(8.5, 5))
    
    names = ['h_promo_count_hzone_cume', 'h_promo_span_hzone_cume', 'h_promo_promo_hzone_cume']
    names_legend = ['Obs. Count', 'Obs. Span', 'Promotion']
    n_plot = len(names)
    
    for n, name in enumerate(names):
        f_mean = f'{name}_q50'
        f_bar1 = f'{name}_q25'
        f_bar2 = f'{name}_q75'
        
        yerr_lower = t_promote[f_mean].values - t_promote[f_bar1].values
        yerr_upper = t_promote[f_bar2].values - t_promote[f_mean].values
        
        ax.errorbar(tsamp + t_offsets[n],
                   t_promote[f_mean].values,
                   yerr=[yerr_lower, yerr_upper],
                   color=line_colors[n],
                   **ebar_props)
    
    style_promote_plot(ax,
                      'Habitable Zone: Cumulative Promotions vs. Detector Time',
                      'Targets Passing Criterion [count]',
                      names_legend)
    write_plots(fig, 'promote-hzone-cume')
    plt.close(fig)
    
    # C: Cumulative Earth plot
    ####################################################################
    
    fig, ax = plt.subplots(figsize=(8.5, 5))
    
    names = ['h_promo_count_earth_cume', 'h_promo_span_earth_cume', 'h_promo_promo_earth_cume']
    names_legend = ['Obs. Count', 'Obs. Span', 'Promotion']
    n_plot = len(names)
    
    for n, name in enumerate(names):
        f_mean = f'{name}_q50'
        f_bar1 = f'{name}_q25'
        f_bar2 = f'{name}_q75'
        
        yerr_lower = t_promote[f_mean].values - t_promote[f_bar1].values
        yerr_upper = t_promote[f_bar2].values - t_promote[f_mean].values
        
        ax.errorbar(tsamp + t_offsets[n],
                   t_promote[f_mean].values,
                   yerr=[yerr_lower, yerr_upper],
                   color=line_colors[n],
                   **ebar_props)
    
    style_promote_plot(ax,
                      'Earthlike Planets: Cumulative Promotions vs. Detector Time',
                      'Targets Passing Criterion [count]',
                      names_legend)
    write_plots(fig, 'promote-earth-cume')
    plt.close(fig)
    
    # D: Cumulative Star plot
    ####################################################################
    
    fig, ax = plt.subplots(figsize=(8.5, 5))
    
    names = ['h_promo_count_star_cume', 'h_promo_span_star_cume', 'h_promo_promo_star_cume']
    names_legend = ['Obs. Count', 'Obs. Span', 'Promotion']
    n_plot = len(names)
    
    for n, name in enumerate(names):
        f_mean = f'{name}_q50'
        f_bar1 = f'{name}_q25'
        f_bar2 = f'{name}_q75'
        
        yerr_lower = t_promote[f_mean].values - t_promote[f_bar1].values
        yerr_upper = t_promote[f_bar2].values - t_promote[f_mean].values
        
        ax.errorbar(tsamp + t_offsets[n],
                   t_promote[f_mean].values,
                   yerr=[yerr_lower, yerr_upper],
                   color=line_colors[n],
                   **ebar_props)
    
    style_promote_plot(ax,
                      'Counted By Star: Earthlike Planets: Cumulative Promotions vs. Detector Time',
                      'Targets Passing Criterion [count]',
                      names_legend)
    write_plots(fig, 'promote-star-cume')
    plt.close(fig)
    
    # E: Cumulative Span options plot
    ####################################################################
    
    fig, ax = plt.subplots(figsize=(8.5, 5))
    
    names = ['h_promo_spanPlan_star_cume', 'h_promo_spanHZ0_star_cume', 
             'h_promo_spanHZ1_star_cume', 'h_promo_spanEarth_star_cume']
    names_legend = ['Span: Any Planet Period', 'Span: Inner HZ Period', 
                    'Span: Outer HZ Period', 'Span: Earthlike Period']
    n_plot = len(names)
    
    for n, name in enumerate(names):
        f_mean = f'{name}_q50'
        f_bar1 = f'{name}_q25'
        f_bar2 = f'{name}_q75'
        
        yerr_lower = t_promote[f_mean].values - t_promote[f_bar1].values
        yerr_upper = t_promote[f_bar2].values - t_promote[f_mean].values
        
        ax.errorbar(tsamp + t_offsets[n],
                    t_promote[f_mean].values,
                    yerr=[yerr_lower, yerr_upper],
                    # color=?
                    **ebar_props)
    
    style_promote_plot(ax,
                      'Counted By Star: Cumulative Observation Span vs. Detector Time',
                      'Targets Passing Span Criterion [count]',
                      names_legend)
    write_plots(fig, 'promote-star-span-cume')
    plt.close(fig)
    
    ####################################################################
    # Promotion-count histograms (function of planet count)
    ####################################################################
    
    # Bail out of all plots if the canary field is not present
    canary = 'h_phist_t1_count_hzone_mean'
    if t_phist.empty or canary not in t_phist.columns:
        print(f'\t{PROGNAME}: No {canary} in promotion table ("phist"), skipping')
        return tracker.get_files()
    
    x_values = t_phist['h_phist_count_lo'].values
    names_legend = ['Obs. Count', 'Obs. Span', 'Promotion']
    
    # Line-styles to distinguish often-overlapping lines
    ls_ct = 'o-'  # ct = count
    ls_sp = 'o-'  # sp = span
    ls_pr = 'x-'  # pr = promo
    ls_x = '+-'   # extra
    
    # Small offset to one line
    del_offset = 0.00
    
    # A: HZ Plot
    ####################################################################
    
    title_tmpl = '%s Planets: Promotion Tiers [%s Mission Time]'
    
    fig, ax = plt.subplots(figsize=(8.5, 5))
    ax.plot(x_values, t_phist['h_phist_t1_count_hzone_mean'].values + del_offset, ls_ct,
            x_values, t_phist['h_phist_t1_span_hzone_mean'].values, ls_sp,
            x_values, t_phist['h_phist_t1_promo_hzone_mean'].values, ls_pr)
    style_phist_plot(ax, title_tmpl % ('Habitable Zone', '2 Years'), names_legend)
    write_plots(fig, 'phist-hzone-2year')
    plt.close(fig)
    
    fig, ax = plt.subplots(figsize=(8.5, 5))
    ax.plot(x_values, t_phist['h_phist_t2_count_hzone_mean'].values + del_offset, ls_ct,
            x_values, t_phist['h_phist_t2_span_hzone_mean'].values, ls_sp,
            x_values, t_phist['h_phist_t2_promo_hzone_mean'].values, ls_pr)
    style_phist_plot(ax, title_tmpl % ('Habitable Zone', '3 Years'), names_legend)
    write_plots(fig, 'phist-hzone-3year')
    plt.close(fig)
    
    # B: Earth plot
    ####################################################################
    
    fig, ax = plt.subplots(figsize=(8.5, 5))
    ax.plot(x_values, t_phist['h_phist_t1_count_earth_mean'].values + del_offset, ls_ct,
            x_values, t_phist['h_phist_t1_span_earth_mean'].values, ls_sp,
            x_values, t_phist['h_phist_t1_promo_earth_mean'].values, ls_pr)
    style_phist_plot(ax, title_tmpl % ('Earthlike', '2 Years'), names_legend)
    write_plots(fig, 'phist-earth-2year')
    plt.close(fig)
    
    fig, ax = plt.subplots(figsize=(8.5, 5))
    ax.plot(x_values, t_phist['h_phist_t2_count_earth_mean'].values + del_offset, ls_ct,
            x_values, t_phist['h_phist_t2_span_earth_mean'].values, ls_sp,
            x_values, t_phist['h_phist_t2_promo_earth_mean'].values, ls_pr)
    style_phist_plot(ax, title_tmpl % ('Earthlike', '3 Years'), names_legend)
    write_plots(fig, 'phist-earth-3year')
    plt.close(fig)
    
    # C: Counting-by-Star plot
    ####################################################################
    
    fig, ax = plt.subplots(figsize=(8.5, 5))
    ax.plot(x_values, t_phist['h_phist_t1_count_star_mean'].values + del_offset, ls_ct,
            x_values, t_phist['h_phist_t1_span_star_mean'].values, ls_sp,
            x_values, t_phist['h_phist_t1_promo_star_mean'].values, ls_pr)
    style_phist_plot(ax, title_tmpl % ('By Star, with Earthlike', '2 Years'), names_legend)
    write_plots(fig, 'phist-star-2year')
    plt.close(fig)
    
    fig, ax = plt.subplots(figsize=(8.5, 5))
    ax.plot(x_values, t_phist['h_phist_t2_count_star_mean'].values + del_offset, ls_ct,
            x_values, t_phist['h_phist_t2_span_star_mean'].values, ls_sp,
            x_values, t_phist['h_phist_t2_promo_star_mean'].values, ls_pr)
    style_phist_plot(ax, title_tmpl % ('By Star, with Earthlike', '3 Years'), names_legend)
    write_plots(fig, 'phist-star-3year')
    plt.close(fig)
    
    # D: Span plot
    ####################################################################
    
    title_tmpl_span = '%s Planets: Various Span Critera Only [%s Mission Time]'
    names_legend_span = ['Obs. Span > T/2 (Any Actual Planet)', 
                        'Obs. Span > T/2 (Inner HZ, Hypothet.)',
                        'Obs. Span > T/2 (Outer HZ, Hypothet.)', 
                        'Obs. Span > T/2 (Actual Earthlike)']
    
    fig, ax = plt.subplots(figsize=(8.5, 5))
    ax.plot(x_values, t_phist['h_phist_t1_spanPlan_star_mean'].values + del_offset, ls_ct,
            x_values, t_phist['h_phist_t1_spanHZ0_star_mean'].values, ls_sp,
            x_values, t_phist['h_phist_t1_spanHZ1_star_mean'].values, ls_pr,
            x_values, t_phist['h_phist_t1_spanEarth_star_mean'].values, ls_x)
    style_phist_plot(ax, title_tmpl_span % ('By Star, with Earthlike', '2 Years'), names_legend_span)
    write_plots(fig, 'phist-star-span-2year')
    plt.close(fig)
    
    fig, ax = plt.subplots(figsize=(8.5, 5))
    ax.plot(x_values, t_phist['h_phist_t2_spanPlan_star_mean'].values + del_offset, ls_ct,
            x_values, t_phist['h_phist_t2_spanHZ0_star_mean'].values, ls_sp,
            x_values, t_phist['h_phist_t2_spanHZ1_star_mean'].values, ls_pr,
            x_values, t_phist['h_phist_t2_spanEarth_star_mean'].values, ls_x)
    style_phist_plot(ax, title_tmpl_span % ('By Star, with Earthlike', '3 Years'), names_legend_span)
    write_plots(fig, 'phist-star-span-3year')
    plt.close(fig)

    return tracker.get_files()


def main():
    """
    Command-line interface for plot_drm_promote
    """
    parser = argparse.ArgumentParser(
        description='Plot target promotion counts in a DRM-set',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
    python plot_drm_promote.py "data/%s.%s" "output/det-%s.%s"
    
The first argument is the source template with two %%s placeholders that will be
filled with ("info", "csv"), ("promote", "csv"), and ("phist", "csv") to 
locate input files.

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
    plot_data = cs.load_csv_files(args.src_tmpl, ['promote', 'promote-hist'])
    rv = plot_drm_promote(reduce_info, plot_data, args.dest_tmpl, mode)
    return rv


if __name__ == '__main__':
    rv = main()
    if rv is None:
        print(f"Plots failed. Error signaled.", file=sys.stderr)
    else:
        print(f"Done. Wrote {len(rv)} plot(s).")
    sys.exit(1 if rv is None else 0)
