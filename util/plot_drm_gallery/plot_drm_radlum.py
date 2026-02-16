#!/usr/bin/env python
"""
Plot radius/luminosity summary for a drm-set

Radius/luminosity plots of detections.
If dest_tmpl is empty, no files are written.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import argparse
import sys
import os
# uniform plot appearance
from . import common_style as cs

# Program name for error messages
PROGNAME = os.path.basename(sys.argv[0])

# Verbosity (also set from mode)
VERBOSE = 1

# overlay imports
import matplotlib.image as mpimg
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from pathlib import Path


def plot_drm_planet_overlay(ax_dest, mode=None):
    """
    Plot a decorative underlay onto ax_dest
    
    Plot a decorative underlay onto ax_dest (typically current axes),
    returning the axis object(s) of the underlay.
    The mode is currently unused.
    
    Parameters
    ----------
    ax_dest : matplotlib.axes.Axes
        Destination axes object where planet overlays will be added
    mode : dict, optional
        Mode dictionary (currently unused)
        
    Returns
    -------
    list of matplotlib.axes.Axes
        List of axis objects for each planet overlay
    """
    
    if mode is None:
        mode = {}
    
    # Configuration selection
    use_config = 3  # 1, 2, or 3
    
    if use_config == 1:
        # 3-class setup (late 2017, early 2018)
        # planets: earthlike, neptune, jupiter, left-to-right
        pfiles = ['gliese-667C_c.png', 'neptune.png', 'jupiter.png']
        # positions are x0, y0, xsize, ysize
        ppos = [[0.15, 0.6, 0.15, 0.15], [0.35, 0.6, 0.3, 0.3], [0.58, 0.45, 0.45, 0.45]]
        # pixels to crop each planet (above) on each edge
        pcrops = [[0, 0, 0, 0], [70, 70, 70, 70], [300, 80, 80, 400]]
    
    elif use_config == 2:
        # 5-class setup (mid 2018)
        # planets, left-to-right: 
        # rocky sub-earth, super-earth, sub-neptune, neptune, jupiter
        pfiles = ['trappist-1_d.png', 'gliese-667C_c.png', 'gj-1214_b.png', 
                  'neptune.png', 'jupiter.png']
        # positions are x0, y0, xsize, ysize
        ppos = [[0.22, 0.4, 0.10, 0.10], [0.33, 0.50, 0.15, 0.15], 
                [0.45, 0.6, 0.2, 0.2], [0.53, 0.62, 0.3, 0.3], [0.63, 0.57, 0.35, 0.35]]
        # pixels to crop each planet on each edge (top, left, bottom, right)
        pcrops = [[25, 25, 25, 25], [0, 0, 0, 0], [0, 0, 0, 0], 
                  [70, 70, 70, 70], [300, 80, 80, 300]]
    
    else:  # use_config == 3
        # 5-class setup, planets to scale (july 2018)
        # planets, left-to-right: 
        # rocky sub-earth, super-earth, sub-neptune, neptune, jupiter (NB: cropped tight)
        pfiles = ['trappist-1_d.png', 'gliese-667C_c.png', 'gj-1214_b.png', 
                  'neptune.png', 'jupiter-tight.png']
        # pixels to crop each planet (above) on each edge (top, left, bottom, right)
        # Note: All crops except Jupiter are to the planetary disk, so that
        # scaling is understood
        pcrops = [[25, 25, 25, 25], [25, 25, 25, 25], [10, 10, 10, 10], 
                  [100, 100, 100, 100], [300, 0, 0, 300]]
        
        j_pre = 0.5  # amount the Jupiter-class image is
                     # already shrunken due to cropping - others are
                     # cropped on the planetary disk
        
        # Set planet positions within the enclosing axis
        # positions are x0, y0, xsize, ysize
        # size is relative to radius property of the corresponding planet class
        ctrline = 0.83  # centerline of the left-to-right line of planets
        size_base = 0.06
        size_rel = np.array([3.0, 4.5, 7.0, 12.0, 22.0*j_pre]) # relative size of planets
        size_img = size_base * size_rel
        ppos = [
            [0.15, ctrline - size_img[0]/2, size_img[0], size_img[0]],
            [0.33, ctrline - size_img[1]/2, size_img[1], size_img[1]],
            [0.41, ctrline - size_img[2]/2, size_img[2], size_img[2]],
            [0.47, ctrline - size_img[3]/2, size_img[3], size_img[3]],
            [0.66, ctrline - size_img[4]/2 + 0.05, size_img[4], size_img[4]]
        ]
    
    # Make planet underlay
    ax_list = []
    
    # Kill white background on current axes, so underlay is visible
    ax_dest.patch.set_facecolor('none')
    
    for p, pfile_ in enumerate(pfiles):
        # FIXME: this temporary location for resources should become canonical
        pfile = 'Local/Matlab/mfile' / Path(pfile_)
        # Check if file exists
        if not os.path.exists(pfile):
            print(f"Warning: Planet image file '{pfile}' not found, skipping")
            continue
        
        # FIXME: this big try/catch is too coarse-grained
        try:
            # Load underlay image
            im_orig = mpimg.imread(pfile)
            
            # Handle alpha channel
            if im_orig.shape[2] == 4:
                # Image has alpha channel
                alfa_orig = im_orig[:, :, 3]
                im_rgb = im_orig[:, :, :3]
            else:
                # No alpha channel, create opaque alpha
                alfa_orig = np.ones(im_orig.shape[:2])
                im_rgb = im_orig
            
            # Crop image and alpha map
            pcrop = pcrops[p]  # crop info, for convenience
            # pcrop = [top, left, bottom, right]
            if pcrop[0] > 0 or pcrop[1] > 0 or pcrop[2] > 0 or pcrop[3] > 0:
                im = im_rgb[pcrop[0]:im_rgb.shape[0]-pcrop[2], 
                           pcrop[1]:im_rgb.shape[1]-pcrop[3], :]
                alfa = alfa_orig[pcrop[0]:alfa_orig.shape[0]-pcrop[2], 
                                pcrop[1]:alfa_orig.shape[1]-pcrop[3]]
            else:
                im = im_rgb
                alfa = alfa_orig
            
            # Create an inset axis for the underlay
            # Position is [x0, y0, width, height] in figure fraction
            pos = ppos[p]
            # I tried 2 ways, an "inset axis" (original, axis object) and
            # a new axis within in the figure (tinkering, figure object)
            if False:
                ax1 = inset_axes(ax_dest,
                             width=f'{pos[2]*100}%', height=f'{pos[3]*100}%',
                             loc='lower left',
                             bbox_to_anchor=(pos[0], pos[1], pos[2], pos[3]),
                             bbox_transform=ax_dest.transAxes,
                             borderpad=0)
            else:
                # in this case it's a figure not an axis, coords are shrunken
                ax1 = ax_dest.add_axes([pos[0], pos[1]*0.75, pos[2]*0.5, pos[3]*0.5],
                                       zorder=3,
                                       facecolor='blue'
                                       )
            
            ### FIXME: frustrated, current status:
            ### -- zorder still not allowing bars above planet images, why?
            ### -- alpha channel of the images is not working end-to-end
            ###       - "alfa" *is* a float, it does change properties
            ###       - but the thing that shows through is black, not the lower axis

            # Plot the image into the axis (including transparency info)
            ax1.imshow(im, alpha=alfa, aspect='equal', zorder=1)
            
            # Set axis parameters
            ax1.set_facecolor('none')  # make the axis itself transparent
            ax1.set_aspect('equal')    # 1:1 aspect ratio so planet-images are circles
            
            # Turn off axis adornments
            ax1.set_xticks([])
            ax1.set_yticks([])
            ax1.axis('off')
            
            # zorder: above rectangles, below plotted bars/lines
            # note double-set with add_axes above
            ax1.set_zorder(5)
            
            # Record the underlay axis, ax1
            ax_list.append(ax1)
            
        except Exception as e:
            print(f"Warning: Could not load planet image '{pfile}': {e}")
            continue
    
    return ax_list


def plot_drm_radlum(reduce_info, src_tmpl, dest_tmpl, mode):
    """
    Plot radius/luminosity summary for a drm-set
    
    Radius/luminosity plots of detections.
    
    Parameters
    ----------
    src_tmpl : str
        Template string for input file paths with two %s placeholders
        Example: "data/%s.%s" will be filled with ("info", "csv"), 
        ("radlum", "csv"), and ("earth", "csv")
    dest_tmpl : str
        Template string for output file paths with two %s placeholders
    mode : dict
        Dictionary with 'op' key containing operation mode string
        
    Outputs
    -------
    Saves plots to disk with names:
        radlum-population.png
        radlum-det.png, radlum-det-all.png
        radlum-char.png, radlum-char-all.png, radlum-char-strict.png
        radlum-char-snr.png
    """
    
    # Load data using the source template
    # update global verbosity
    global VERBOSE
    VERBOSE = mode.get('verbose', VERBOSE)

    try:
        radlum_file = src_tmpl % ("radlum", "csv")
        t_radlum = pd.read_csv(radlum_file)
    except Exception as e:
        print(f"{PROGNAME}: Fatal: Could not load radlum file: {e}", file=sys.stderr)
        sys.exit(1)
    
    try:
        earth_file = src_tmpl % ("earth", "csv")
        t_earth = pd.read_csv(earth_file)
    except Exception as e:
        print(f"{PROGNAME}: Fatal: Could not load earth file: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Skip unless mode.op contains our name or a *
    if '*' not in mode.get('op', '') and 'radlum' not in mode.get('op', ''):
        print('Radius/luminosity plots: skipping, as directed.')
        return []
    
    ##################################################################
    # Common data and functions
    ##################################################################
    
    # Number of radii obtained by the bin boundaries
    n_rad = len(np.unique(t_radlum['h_Rp_lo'].values))
    # Number of luminosities
    # Obtain from the number of total bins, divided by n_rad
    # Cannot use #uniques because L bin boundaries vary by radius
    n_lum = len(t_radlum['h_L_lo'].values) // n_rad
    
    # X-axis with larger inter-radius spacing
    # e.g.: [[1 2 3] [4.2 5.2 6.2] [7.4 8.4 9.4] [10.6 11.6 12.6]]
    x_bin_xtra = 0.29
    x_bin = np.reshape(
        np.arange(1, n_lum * n_rad + 1).reshape(n_rad, n_lum) + 
        np.arange(n_rad).reshape(n_rad, 1) * x_bin_xtra,
        -1
    )
    x_bin_plus = np.concatenate([[-x_bin_xtra], x_bin])  # includes earths
    
    # Bar color (RGB triple) for each luminosity: hot, warm, cold
    colors = [
        [0.8667, 0.1843, 0.1020],
        [0.3686, 0.6706, 0.9608],
        [0.8275, 0.9412, 1.0]
    ]
    
    # For axis labels
    lumens = ['Hot', 'Warm', 'Cold']
    lumens_plus = ['Earth'] + lumens * n_rad
    
    # Error-bar color and style
    color_eb = [0.3, 0.3, 0.3]
    
    # Text-block style
    # Appears on top, but with high bbox transparency
    style_block = {
        'zorder': 10, # on top
        'alpha': 0.8, # alpha for text only
        'fontsize': 10,
        'verticalalignment': 'top',
        'bbox': dict(boxstyle='round', alpha=0.5, facecolor=[0.8, 0.8, 0.8], edgecolor='k')
    }
    
    # Helper matrix for stacked bar plots
    # Collapses a set of counts down
    # Makes a 9x9 diagonal into a 9x3 set of stacked plots
    collapse = np.tile(np.eye(n_lum), (n_rad, 1))
    
    # File extensions to write
    ext_list = ['png']

    # Track output files
    tracker = cs.PlotTracker()
    tracker.set_ext_list(ext_list)

    # Inner function: rectangle underlay
    def place_rect_underlay(ax, n_rad, n_lum, x_bin, x_bin_xtra):
        """Place vertical stripes to visually group planets of same radius class"""
        # Save plot height and width, for the rectangles below
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        
        # Add Earth (0) plus n_rad (1..n_rad) rectangles
        # - one for each radius class
        for k1 in range(n_rad + 1):
            # Rectangle is between x1 and x2, and top-to-bottom
            if k1 == 0:
                # Special-case Earth
                x1 = x_bin_plus[0] - 0.5 - 0.2  # 0.5 = half-bin; 0.2 = empirical bar width
                x2 = x_bin[0] - 0.5 * (x_bin_xtra + 1)  # end: left of first x_bin
                fc = [0.7569, 0.8863, 0.7765]  # pale green
            else:
                fc = [0.85, 0.85, 0.85]  # pale gray
                # Index into bins: 0, n_lum, ...
                k2 = (k1 - 1) * n_lum
                # We pad by 1 + x_bin_xtra between radius groups
                x1 = x_bin[k2] - 0.5 * (x_bin_xtra + 1)
                x2 = x1 + n_lum + x_bin_xtra
                if x2 > xlim[1]:
                    print(f'\tvanishing rect {k1}') # FIXME/REMOVE
                x2 = min(x2, xlim[1])  # cap rightmost edge at plot end
            
            # Plot the rectangle itself
            x_gap = 0.11  # (half the) gap between adjacent rectangular panels
            rect = patches.Rectangle(
                (x1 + x_gap, ylim[0]), # SW corner
                (x2 - x1) - 2 * x_gap, # width
                ylim[1] - ylim[0], # height
                facecolor=fc,
                linestyle='none',
                zorder=0
            )
            ax.add_patch(rect)
    
    # Inner function: write the current figure to files
    def write_plots(fig, dest_name):
        """Write the current figure to various files"""
        # FIXME: facecolor disabled because of overlay confusion
        tracker.write_plots(fig, dest_name, dest_tmpl, verbose=VERBOSE, facecolor=None)

    # Inner function: Set up plot/axis styles, title, axis labels
    def style_bar_plot(ax, obs_type):
        """Style the bar plot with title, labels, and formatting"""
        ax.grid(True, axis='y')
        ax.set_xticks(x_bin_plus)
        ax.set_xticklabels(lumens_plus)
        ax.set_xlim([x_bin_plus[0] - 1, x_bin[-1] + 1])
        
        # Errorbars can go < 0
        ylim = ax.get_ylim()
        ax.set_ylim(max(0, ylim[0]), ylim[1])
        
        title1 = f'{obs_type} vs. Insolation and Planet Radius'
        title2 = cs.plot_make_title(reduce_info)
        
        ax.set_title(f'{title2}\n{title1}', fontsize=11*1.1, fontweight='bold')
        ax.set_xlabel('Planet Type and Size', fontweight='bold')
        ax.set_ylabel(obs_type, fontweight='bold')
        ax.tick_params(labelsize=10)
    
    # uniform errorbar() properties
    ebar_props = {
        "linewidth": 1.5,
        "markersize": 5,
        "capsize": 2,
        "color": color_eb,
        "zorder": 11 # on top
        }
    # position of all text blocks
    text_pos = (0.01, 0.98)
    # uniform bar() style -- on top
    bar_style = {'width': 0.8, 'zorder': 10}

    ##################################################################
    # Overall Planet Population
    ##################################################################
    
    # Data, for convenience
    h_counts = t_radlum['h_RpL_population_mean'].values
    h_stds = t_radlum['h_RpL_population_std'].values
    
    # Basic plot method: use a bar-plot with overlaid error bars
    # Separates h_counts into n_rad groups which will have repeating colors
    h_grouped = np.diag(h_counts) @ collapse
    
    fig, ax = plt.subplots(figsize=(9.5, 5))
    
    # Basic bar plot - stacked bars
    bottom = np.zeros(len(x_bin))
    for k in range(h_grouped.shape[1]):
        color_idx = k % len(colors)
        ax.bar(x_bin, h_grouped[:, k], bottom=bottom, 
               color=colors[color_idx], **bar_style)
        bottom += h_grouped[:, k]
    
    # Exo-earth barplot and errorbar
    ct_e = t_earth['exoE_population_mean'].values[0]
    eb_e = t_earth['exoE_population_std'].values[0]
    ax.bar([-x_bin_xtra], [ct_e], color=[0.1, 0.7, 0.3], **bar_style)
    ax.errorbar(-x_bin_xtra, ct_e, eb_e, fmt='.', 
                **ebar_props)
    
    # First-stack errorbars
    ax.errorbar(x_bin, h_counts, h_stds, fmt='.', 
                **ebar_props)
    
    # Establish plot/axis styles + labels
    style_bar_plot(ax, 'Planet Occurrence Rate')
    
    # Explanatory legend
    block_text = ('Counting Over All Planets in Target List\n'
                 'Normalization: Planets of that Type, per Star\n'
                 'Earth Column is Observed $\eta$ Earth\n'
                 'Error Bar: $\pm$ 1 sigma')
    ax.text(*text_pos, block_text, transform=ax.transAxes, **style_block)
    
    # Planet overlay
    # (Disabled)
    #   Original call from claude:
    #    ax_ulay = plot_drm_planet_overlay(ax, mode)
    #   Tinkering with an alternate axis setup:
    #    ax_ulay = plot_drm_planet_overlay(fig, mode)
    #
    # Rectangle underlay
    place_rect_underlay(ax, n_rad, n_lum, x_bin, x_bin_xtra)
    
    # Write the plots
    write_plots(fig, 'radlum-population')
    plt.close(fig)
    
    ##################################################################
    # Detections - unique
    ##################################################################
    
    # Data, for convenience
    h_counts = t_radlum['h_RpL_det_main_mean'].values
    h_stds = t_radlum['h_RpL_det_main_std'].values
    h_counts_alt = t_radlum['h_RpL_det_alt_mean'].values
    h_stds_alt = t_radlum['h_RpL_det_alt_std'].values
    # Pooled standard error of the sum: std = sqrt(std^2 + std^2)
    h_stds_alt2 = np.hypot(h_stds, h_stds_alt)
    
    # Basic plot method
    h_grouped = np.diag(h_counts) @ collapse
    h_grouped_alt = np.diag(h_counts_alt) @ collapse
    
    fig, ax = plt.subplots(figsize=(9.5, 5))
    
    # Basic bar plot - stacked bars (both main and alt)
    bottom = np.zeros(len(x_bin))
    for k in range(h_grouped.shape[1]):
        color_idx = k % len(colors)
        ax.bar(x_bin, h_grouped[:, k], bottom=bottom, 
               color=colors[color_idx], **bar_style)
        bottom += h_grouped[:, k]
    for k in range(h_grouped_alt.shape[1]):
        color_idx = k % len(colors)
        ax.bar(x_bin, h_grouped_alt[:, k], bottom=bottom, 
               color=colors[color_idx], **bar_style)
        bottom += h_grouped_alt[:, k]
    
    # Exo-earth barplot and errorbar
    ct_e = np.array([t_earth['exoE_det_main_mean'].values[0], 
                     t_earth['exoE_det_alt_mean'].values[0]])
    eb_e = np.hypot(t_earth['exoE_det_main_std'].values[0], 
                    t_earth['exoE_det_alt_std'].values[0])
    # Stacked bars for earth
    ax.bar([-x_bin_xtra], [ct_e[0]], color=[0.1, 0.7, 0.3], **bar_style)
    ax.bar([-x_bin_xtra], [ct_e[1]], bottom=[ct_e[0]], 
           color=[0.1, 0.7, 0.3], **bar_style)
    ax.errorbar(-x_bin_xtra, np.sum(ct_e), eb_e, fmt='.', 
                **ebar_props)
    
    # First-stack errorbars
    ax.errorbar(x_bin, h_counts, h_stds, fmt='.', 
                **ebar_props)
    # Second-stack errorbars
    ax.errorbar(x_bin - 0.2, h_counts + h_counts_alt, h_stds_alt, fmt='.', 
                **ebar_props)
    # Pooled errorbars
    ax.errorbar(x_bin + 0.2, h_counts + h_counts_alt, h_stds_alt2, fmt='.', 
                **ebar_props)
    
    # Establish plot/axis styles + labels
    style_bar_plot(ax, 'Unique Detections')
    
    # Explanatory legend
    block_text = ('Lower Segment: Blue Detection Mode\n'
                 'Upper Segment: Combo Detection Mode\n'
                 'Error Bar: $\pm$ 1 sigma')
    ax.text(*text_pos, block_text, transform=ax.transAxes, **style_block)
    
    # Planet + rectangle underlays
    place_rect_underlay(ax, n_rad, n_lum, x_bin, x_bin_xtra)
    
    # Write the plots
    write_plots(fig, 'radlum-det')
    plt.close(fig)
    
    ##################################################################
    # Detections - all
    ##################################################################
    
    # Data, for convenience
    h_counts = t_radlum['h_RpL_xdet_main_mean'].values
    h_stds = t_radlum['h_RpL_xdet_main_std'].values
    h_counts_alt = t_radlum['h_RpL_xdet_alt_mean'].values
    h_stds_alt = t_radlum['h_RpL_xdet_alt_std'].values
    h_stds_alt2 = np.hypot(h_stds, h_stds_alt)
    
    h_grouped = np.diag(h_counts) @ collapse
    h_grouped_alt = np.diag(h_counts_alt) @ collapse
    
    fig, ax = plt.subplots(figsize=(9.5, 5))
    
    # Basic bar plot
    bottom = np.zeros(len(x_bin))
    for k in range(h_grouped.shape[1]):
        color_idx = k % len(colors)
        ax.bar(x_bin, h_grouped[:, k], bottom=bottom, 
               color=colors[color_idx], **bar_style)
        bottom += h_grouped[:, k]
    for k in range(h_grouped_alt.shape[1]):
        color_idx = k % len(colors)
        ax.bar(x_bin, h_grouped_alt[:, k], bottom=bottom, 
               color=colors[color_idx], **bar_style)
        bottom += h_grouped_alt[:, k]
    
    # Exo-earth barplot and errorbar
    ct_e = np.array([t_earth['exoE_xdet_main_mean'].values[0], 
                     t_earth['exoE_xdet_alt_mean'].values[0]])
    eb_e = np.hypot(t_earth['exoE_xdet_main_std'].values[0], 
                    t_earth['exoE_xdet_alt_std'].values[0])
    ax.bar([-x_bin_xtra], [ct_e[0]], color=[0.1, 0.7, 0.3], **bar_style)
    ax.bar([-x_bin_xtra], [ct_e[1]], bottom=[ct_e[0]], 
           color=[0.1, 0.7, 0.3], **bar_style)
    ax.errorbar(-x_bin_xtra, np.sum(ct_e), eb_e, fmt='.', 
                **ebar_props)
    
    # Errorbars
    ax.errorbar(x_bin, h_counts, h_stds, fmt='.', 
                **ebar_props)
    ax.errorbar(x_bin - 0.2, h_counts + h_counts_alt, h_stds_alt, fmt='.', 
                **ebar_props)
    ax.errorbar(x_bin + 0.2, h_counts + h_counts_alt, h_stds_alt2, fmt='.', 
                **ebar_props)
    
    # Establish plot/axis styles + labels
    style_bar_plot(ax, 'All Detections [With Repeats]')
    
    # Explanatory legend
    block_text = ('Lower Segment: Blue Detection Mode\n'
                 'Upper Segment: Combo Detection Mode\n'
                 'Error Bar: $\pm$ 1 sigma')
    ax.text(*text_pos, block_text, transform=ax.transAxes, **style_block)
    
    # Planet + rectangle underlays
    place_rect_underlay(ax, n_rad, n_lum, x_bin, x_bin_xtra)
    
    # Write the plots
    write_plots(fig, 'radlum-det-all')
    plt.close(fig)
    
    ##################################################################
    # Characterizations
    ##################################################################
    
    # Controls the characterization plots we make
    char_plot_menu = [
        ['', '_', '', 'dual', 'radlum-char'],
        ['x', '_', 'With Repeats', 'dual', 'radlum-char-all'],
        ['', 'strict_', 'Strict', 'solo', 'radlum-char-strict'],
    ]
    
    for plot_num, plot_spec in enumerate(char_plot_menu):
        if len(plot_spec) != 5:
            print(f"{PROGNAME}: Fatal: Plot selector {plot_num} has wrong number of elements",
                  file=sys.stderr)
            sys.exit(1)
        
        plot_prop_x = plot_spec[0]
        plot_prop = plot_spec[1]
        plot_title_raw = plot_spec[2]
        plot_mode = plot_spec[3]
        plot_file = plot_spec[4]
        
        if not plot_title_raw:
            plot_title = ''
        else:
            plot_title = f' [{plot_title_raw}]'
        
        # Not all bands will necessarily be present in the t_radlum table
        if plot_mode == 'dual':
            plot_prop_root = '_full'
        else:
            plot_prop_root = '_'  # (strict)
        
        canary = f'h_RpL_{plot_prop_x}char{plot_prop_root}{plot_prop}mean'
        if canary not in t_radlum.columns:
            print(f'{PROGNAME}: Property {canary} not in CSV, skipping the {plot_file} plot.')
            continue
        
        # Data, for convenience
        h_counts = t_radlum[f'h_RpL_{plot_prop_x}char{plot_prop_root}{plot_prop}mean'].values
        h_stds = t_radlum[f'h_RpL_{plot_prop_x}char{plot_prop_root}{plot_prop}std'].values
        
        if plot_mode == 'dual':
            h_counts_alt = t_radlum[f'h_RpL_{plot_prop_x}char_part{plot_prop}mean'].values
            h_stds_alt = t_radlum[f'h_RpL_{plot_prop_x}char_part{plot_prop}std'].values
            h_stds_alt2 = np.hypot(h_stds, h_stds_alt)
        else:
            h_counts_alt = None
            h_stds_alt = None
            h_stds_alt2 = None
        
        # Basic plot method
        h_grouped = np.diag(h_counts) @ collapse
        if h_counts_alt is not None:
            h_grouped_alt = np.diag(h_counts_alt) @ collapse
        else:
            h_grouped_alt = None
        
        fig, ax = plt.subplots(figsize=(9.5, 5))
        
        # Basic bar plot
        bottom = np.zeros(len(x_bin))
        for k in range(h_grouped.shape[1]):
            color_idx = k % len(colors)
            ax.bar(x_bin, h_grouped[:, k], bottom=bottom, 
                   color=colors[color_idx], **bar_style)
            bottom += h_grouped[:, k]
        
        if h_grouped_alt is not None:
            for k in range(h_grouped_alt.shape[1]):
                color_idx = k % len(colors)
                ax.bar(x_bin, h_grouped_alt[:, k], bottom=bottom, 
                       color=colors[color_idx], **bar_style)
                bottom += h_grouped_alt[:, k]
        
        # Exo-earth barplot and errorbar
        if plot_mode == 'dual':
            ct_e = np.array([
                t_earth[f'exoE_{plot_prop_x}char_full{plot_prop}mean'].values[0],
                t_earth[f'exoE_{plot_prop_x}char_part{plot_prop}mean'].values[0]
            ])
            eb_e = np.hypot(
                t_earth[f'exoE_{plot_prop_x}char_full{plot_prop}std'].values[0],
                t_earth[f'exoE_{plot_prop_x}char_part{plot_prop}std'].values[0]
            )
            ax.bar([-x_bin_xtra], [ct_e[0]], color=[0.1, 0.7, 0.3], **bar_style)
            ax.bar([-x_bin_xtra], [ct_e[1]], bottom=[ct_e[0]], 
                   color=[0.1, 0.7, 0.3], **bar_style)
            ax.errorbar(-x_bin_xtra, np.sum(ct_e), eb_e, fmt='.', 
                        **ebar_props)
        else:
            ct_e = t_earth[f'exoE_{plot_prop_x}char_{plot_prop}mean'].values[0]
            eb_e = t_earth[f'exoE_{plot_prop_x}char_{plot_prop}std'].values[0]
            ax.bar([-x_bin_xtra], [ct_e], color=[0.1, 0.7, 0.3], **bar_style)
            ax.errorbar(-x_bin_xtra, ct_e, eb_e, fmt='.', 
                        **ebar_props)
        
        # First-stack errorbars
        ax.errorbar(x_bin, h_counts, h_stds, fmt='.', 
                    **ebar_props)
        
        if h_counts_alt is not None:
            # Second-stack errorbars
            ax.errorbar(x_bin - 0.2, h_counts + h_counts_alt, h_stds_alt, fmt='.', 
                        **ebar_props)
            # Combo errorbars
            ax.errorbar(x_bin + 0.2, h_counts + h_counts_alt, h_stds_alt2, fmt='.', 
                        **ebar_props)
        
        # Establish plot/axis styles + labels
        qualifier = 'All' if plot_prop_x else 'Unique'
        style_bar_plot(ax, f'{qualifier} Characterizations{plot_title}')
        
        # Explanatory legend
        block_top = f'{qualifier} Characterizations{plot_title}'
        block_eb = '  Error Bar: $\pm$ 1 sigma'
        if h_counts_alt is not None:
            block_text = (f'{block_top}\n'
                         '  Lower Segment: Full Characterization\n'
                         '  Upper Segment: Partial Characterization\n'
                         f'{block_eb}')
        else:
            block_text = f'{block_top}\n{block_eb}'
        
        ax.text(*text_pos, block_text, transform=ax.transAxes, **style_block)
        
        # Planet + rectangle underlays
        place_rect_underlay(ax, n_rad, n_lum, x_bin, x_bin_xtra)
        
        # Write the plots
        write_plots(fig, plot_file)
        plt.close(fig)
    
    ##################################################################
    # Characterization SNR
    ##################################################################
    
    snr_plot_menu = [
        ['_', '', 'radlum-char-snr'],
    ]
    # marker adjustment because SNR plot is different
    ebar_props_snr = dict(ebar_props, markersize=7, markerfacecolor='tab:blue')

    for plot_num, plot_spec in enumerate(snr_plot_menu):
        if len(plot_spec) != 3:
            print(f"{PROGNAME}: Fatal: Plot selector {plot_num} has wrong number of elements",
                  file=sys.stderr)
            sys.exit(1)
        
        plot_prop = plot_spec[0]
        plot_title_raw = plot_spec[1]
        plot_file = plot_spec[2]
        
        if not plot_title_raw:
            plot_title = ''
        else:
            plot_title = f' [{plot_title_raw}]'
        
        # Not all bands will necessarily be present in the table
        test_name = f'h_RpL_char_full{plot_prop}mean'
        if test_name not in t_radlum.columns:
            continue
        
        fig, ax = plt.subplots(figsize=(9.5, 5))
        
        # Characterization SNR plot
        # Quantile error bar (asymmetric)
        q50 = t_radlum[f'h_RpL_char_snr{plot_prop}q50'].values
        q25 = t_radlum[f'h_RpL_char_snr{plot_prop}q25'].values
        q75 = t_radlum[f'h_RpL_char_snr{plot_prop}q75'].values
        
        ax.errorbar(x_bin, q50,
                   yerr=[q50 - q25, q75 - q50],
                   fmt='o',
                   **ebar_props_snr)
        
        # Add Earth
        q50_e = t_earth[f'exoE_char_snr{plot_prop}q50'].values[0]
        q25_e = t_earth[f'exoE_char_snr{plot_prop}q25'].values[0]
        q75_e = t_earth[f'exoE_char_snr{plot_prop}q75'].values[0]
        
        ax.errorbar(-x_bin_xtra, q50_e,
                    yerr=[[q50_e - q25_e], [q75_e - q50_e]],
                    fmt='o',
                    # markersize=10, 
                    **ebar_props)
        
        # Limit plot y-range
        # (A) SNR cannot be negative, but error bars sometimes go < 0
        # (B) SNR upper YLim can sometimes be much too high
        top_snr = 1.1 * np.max(t_radlum[f'h_RpL_char_snr{plot_prop}mean'].values)
        ylim_snr = ax.get_ylim()
        ax.set_ylim(max(0, ylim_snr[0]), min(ylim_snr[1], top_snr))
        
        # Establish plot/axis styles + labels
        style_bar_plot(ax, f'Characterization SNR{plot_title}')
        
        # Explanatory legend
        block_text = (f'Characterization SNR{plot_title}\n'
                     '  Circles: Median\n'
                     '  Error Bars: P25, P75')
        ax.text(*text_pos, block_text, transform=ax.transAxes, **style_block)
        
        # Planet + rectangle underlays
        place_rect_underlay(ax, n_rad, n_lum, x_bin, x_bin_xtra)
        
        # Put the axis grid on top of the above annotations
        ax.set_axisbelow(False)
        
        # Write the plots
        write_plots(fig, plot_file)
        plt.close(fig)

    return tracker.get_files()


def main():
    """
    Command-line interface for plot_drm_radlum
    """
    parser = argparse.ArgumentParser(
        description='Plot radius/luminosity summary for a DRM-set',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
    python plot_drm_radlum.py "data/%s.%s" "output/det-%s.%s"
    
The first argument is the source template with two %%s placeholders that will be
filled with ("info", "csv"), ("radlum", "csv"), and ("earth", "csv") to 
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
        plot_drm_radlum(reduce_info, args.src_tmpl, args.dest_tmpl, mode)
    except Exception as e:
        print(f"{PROGNAME}: Fatal: Unexpected error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
