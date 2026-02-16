#!/usr/bin/env python
"""
Plot attempted earth characterizations in a drm-set

Plots of attempted characterizations, e.g., fails/successes in
WA/dMag coordinates.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
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

# Verbosity (also set from mode)
VERBOSE = 1

def plot_drm_earth_chars(reduce_info, plot_data, dest_tmpl, mode):
    """
    Plot attempted earth characterizations in a drm-set

    Plots of attempted characterizations, e.g., fails/successes in
    WA/dMag coordinates.

    Parameters
    ----------
    reduce_info : dict
        Metadata dict from reduce-info.csv
    plot_data : list of DataFrame
        Pre-loaded CSV data [earth-char-list]
    dest_tmpl : str
        Template string for output file paths with two %s placeholders
    mode : dict
        Dictionary with 'op' key containing operation mode string

    Outputs
    -------
    Saves plots to disk with names:
        earth-char-wa-dmag-promo-snr.png
        earth-char-wa-dmag-promo-phi.png
        earth-char-wa-dmag-promo-vmag.png
        earth-char-wa-dmag-promo-zzz-dmag.png
        earth-char-hist-phi.png
        earth-char-hist-phi-log.png
    """

    # update global verbosity
    global VERBOSE
    VERBOSE = mode.get('verbose', VERBOSE)

    # Unpack CSV data
    t_earth_chars, = plot_data
    
    # Skip unless mode.op contains our name or a *
    if '*' not in mode.get('op', '') and 'earth-chars' not in mode.get('op', ''):
        print('Earth char attempt plots: skipping, as directed.')
        return []
    
    # File extensions to write
    ext_list = ['png']

    # Track output files
    tracker = cs.PlotTracker()
    tracker.set_ext_list(ext_list)

    # Helper function to create title from t_info
    # def plot_make_title(t_info):
    #     """Create plot title from metadata"""
    #     if t_info.empty:
    #         rv = ''
    #     elif 'experiment' in t_info:
    #         exp_name = str.strip(t_info['experiment'].iloc[0])
    #         if len(exp_name) < 50:
    #             chaser = f", Ensemble Size {t_info['ensemble_size'].iloc[0]}"
    #         else:
    #             chaser = ''
    #         rv = exp_name + chaser
    #     else:
    #         rv = ''
    #     return rv
    
    # Inner function: Set up plot/axis styles, title, axis labels
    def style_wa_dmag_plot(ax, title1, bartext):
        """Style the WA/dMag plot with title, labels, and colorbar"""
        # Clip Y range at 0
        ylim = ax.get_ylim()
        ax.set_ylim(max(0, ylim[0]), ylim[1])
        
        # Format the title
        title2 = cs.plot_make_title(reduce_info)
        
        # Set title (preventing special interpretation of _) with bold
        ax.set_title(f'{title2}\n{title1}', fontsize=11*1.1, fontweight='bold')
        
        # Axis labels
        ax.set_xlabel('Working Angle (mas)', fontweight='bold')
        ax.set_ylabel('Delta Magnitude', fontweight='bold')
        ax.tick_params(labelsize=13)
        ax.grid(True)
        
        # Colorbar
        if len(bartext) > 0:
            cbar = plt.colorbar(ax.collections[-1], ax=ax)
            cbar.set_label(bartext, fontweight='bold')
        
        # "Easy" rectangle:
        # upper-left corner at: (x=70 mas, y=26 mag)
        # lower-left corner at  (x=70 mas, y=ax[2])
        ax_lim = ax.axis()
        rect = patches.Rectangle((70, ax_lim[2]), ax_lim[1] - 70, 26 - ax_lim[2],
                                 linewidth=1, linestyle='--', 
                                 edgecolor='red', facecolor='none')
        ax.add_patch(rect)
    
    # Inner function: write the current figure to files
    def write_plots(fig, dest_name):
        """Write the current figure to various files"""
        tracker.write_plots(fig, dest_name, dest_tmpl, verbose=VERBOSE, dpi=250)

    # For missions without chars, the table will be basically empty.
    # Skip such tables, there is nothing to do.
    if 'is_success' not in t_earth_chars.columns:
        print('\tNo char info, skipping earth chars plots')
        return []
    
    # Booleans
    ok = t_earth_chars['is_success'].values > 0
    ok2 = (t_earth_chars['is_success'].values == 0) & (t_earth_chars['n_success'].values > 0)
    fail = ~(ok | ok2)
    deep = t_earth_chars['is_deep'].values > 0
    promo = (t_earth_chars['is_promo'].values > 0) & (~deep)
    
    # Real numbers
    ec_wa = t_earth_chars['WA'].values
    ec_dmag = t_earth_chars['dMag'].values
    ec_snr = t_earth_chars['char_SNR'].values
    ec_phi = t_earth_chars['phi'].values
    
    # Not all generated CSVs have this field (3/2019)
    if 'MV' in t_earth_chars.columns:
        ec_vmag = t_earth_chars['MV'].values
    else:
        ec_vmag = np.zeros_like(ec_wa)
    
    # Derived
    easy = (ec_dmag < 26) & (ec_wa > 70)
    
    # Table of wa/dmag plots to make.
    # Format: textual full name, abbreviated name, index variable
    dprops = {'s': 300}
    pprops = {'s': 100, 'alpha': 0.5}
    
    # Eliminated deep-dive plots, 2024/12
    plot_roster = [
        ['Promoted', 'promo', promo, pprops]
    ]
    
    # Create new figures for each plot?
    new_figs = True
    
    for pnum, plot_spec in enumerate(plot_roster):
        full_name = plot_spec[0]
        abbr_name = plot_spec[1]
        c_inx = plot_spec[2]
        gprops = plot_spec[3]
        
        ################################################################
        # WA/dMag, shading = char SNR
        ################################################################
        
        fig, ax = plt.subplots(figsize=(8.5, 5))
        
        ax.scatter(ec_wa[ok], ec_dmag[ok], s=10, c=[[0.5, 0.5, 0.5]], marker='+')
        ps = ax.scatter(ec_wa[fail & c_inx], ec_dmag[fail & c_inx], 
                       c=ec_snr[fail & c_inx], marker='.', **gprops)
        ps2 = ax.scatter(ec_wa[ok2 & c_inx], ec_dmag[ok2 & c_inx], s=10,
                        c=ec_snr[ok2 & c_inx], marker='x')
        
        ax_lim = ax.axis()
        txt_x = (ax_lim[0] + ax_lim[1]) / 2  # x midline for text annotations
        txt_y = 0.7 * ax_lim[2] + 0.3 * ax_lim[3]  # y midline for text annotations
        
        style_wa_dmag_plot(ax,
            f'{full_name}: Earth Characterizations vs. WA and dMag, shaded by Char. SNR',
            'Char SNR')
        ax.legend([f'Successful Chars ({np.sum(ok)})',
                  f'Failed {full_name} Chars ({np.sum(fail & c_inx)})',
                  f'Split Multi-Earth {full_name} Chars ({np.sum(ok2 & c_inx)})'],
                 loc='upper right')
        
        # Easy chars sub-legend
        ax.text(txt_x, txt_y,
               f'"Easy" Characterizations within red box:\n'
               f'  {np.sum(easy & ok)} Successful\n'
               f'  {np.sum(easy & c_inx & fail)} Failed (counting {full_name} only)',
               fontsize=13)
        
        write_plots(fig, f'earth-char-wa-dmag-{abbr_name}-snr')
        plt.close(fig)
        
        ################################################################
        # WA/dMag, shading = log10(Phi)
        ################################################################
        
        fig, ax = plt.subplots(figsize=(8.5, 5))
        
        ax.scatter(ec_wa[ok], ec_dmag[ok], s=10, c=[[0.5, 0.5, 0.5]], marker='+')
        ps = ax.scatter(ec_wa[~ok & c_inx], ec_dmag[~ok & c_inx],
                       c=np.log10(ec_phi[~ok & c_inx]), marker='.', **gprops)
        
        ax_lim = ax.axis()
        
        style_wa_dmag_plot(ax,
            f'{full_name}: Earth Characterizations vs. WA and dMag, shaded by Log Phi',
            'log$_{{10}}$(?) : truncated at -2')
        ax.legend([f'Successful Chars ({np.sum(ok)})',
                  f'Failed {full_name} Chars ({np.sum(~ok & c_inx)})'],
                 loc='upper right')
        
        # Easy chars sub-legend
        txt_x = (ax_lim[0] + ax_lim[1]) / 2
        txt_y = 0.7 * ax_lim[2] + 0.3 * ax_lim[3]
        ax.text(txt_x, txt_y,
               f'"Easy" Characterizations within red box:\n'
               f'  {np.sum(easy & ok)} Successful\n'
               f'  {np.sum(easy & c_inx & ~ok)} Failed (counting {full_name} only)',
               fontsize=13)
        
        # -2 on log10 scale = 0.01 = 5 magnitudes
        ps.set_clim(-2, 0)
        
        write_plots(fig, f'earth-char-wa-dmag-{abbr_name}-phi')
        plt.close(fig)
        
        ################################################################
        # WA/dMag, shading = Vmag
        ################################################################
        
        fig, ax = plt.subplots(figsize=(8.5, 5))
        
        ax.scatter(ec_wa[ok], ec_dmag[ok], s=10, c=[[0.5, 0.5, 0.5]], marker='+')
        ps = ax.scatter(ec_wa[~ok & c_inx], ec_dmag[~ok & c_inx],
                       c=ec_vmag[~ok & c_inx] + ec_dmag[~ok & c_inx], 
                       marker='.', **gprops)
        
        ax_lim = ax.axis()
        
        style_wa_dmag_plot(ax,
            f'{full_name}: Earth Characterizations vs. WA and dMag, shaded by Planet Magnitude',
            'Planet Apparent Visual Magnitude')
        ax.legend([f'Successful Chars ({np.sum(ok)})',
                  f'Failed {full_name} Chars ({np.sum(~ok & c_inx)})'],
                 loc='upper right')
        
        # Easy chars sub-legend
        txt_x = (ax_lim[0] + ax_lim[1]) / 2
        txt_y = 0.7 * ax_lim[2] + 0.3 * ax_lim[3]
        ax.text(txt_x, txt_y,
               f'"Easy" Characterizations within red box:\n'
               f'  {np.sum(easy & ok)} Successful\n'
               f'  {np.sum(easy & c_inx & ~ok)} Failed (counting {full_name} only)',
               fontsize=13)
        
        write_plots(fig, f'earth-char-wa-dmag-{abbr_name}-vmag')
        plt.close(fig)
        
        ################################################################
        # WA/dMag, shading = Adjusted Dmag
        ################################################################
        
        fig, ax = plt.subplots(figsize=(8.5, 5))
        
        # Adjusted dmag - improvement possible
        ec_dmagPhi = ec_dmag + np.maximum(-5, 2.5 * np.log10(np.minimum(0.7, ec_phi) / 0.7))
        
        ax.scatter(ec_wa[ok], ec_dmag[ok], s=10, c=[[0.5, 0.5, 0.5]], marker='+')
        ps = ax.scatter(ec_wa[~ok & c_inx], ec_dmagPhi[~ok & c_inx],
                       c=ec_phi[~ok & c_inx], marker='.', **gprops)
        
        style_wa_dmag_plot(ax,
            f'{full_name}: Earth Characterizations vs. WA and ADJUSTED dMag, shaded by Phi',
            'Lambertian ?')
        ax.legend(['Successful Chars', f'Failed {full_name} Chars'],
                 loc='upper right')
        
        write_plots(fig, f'earth-char-wa-dmag-{abbr_name}-zzz-dmag')
        plt.close(fig)
    
    ####################################################################
    # Remaining plots
    ####################################################################
    
    # Supplementary title
    title_x = cs.plot_make_title(reduce_info)
    
    tprops = {'fontweight': 'bold'}
    
    ################################################################
    # Histogram of Phi
    ################################################################
    
    fig = plt.figure(figsize=(8.5, 5))
    
    # Brighten the colormap
    cmap = plt.cm.jet(np.linspace(0, 1, 256))
    cmap = 0.5 + 0.5 * cmap  # brighten
    
    h_bins = np.arange(0, 1.04, 0.04)
    h_promo_fail, _ = np.histogram(ec_phi[~ok & promo], h_bins)
    h_deep_fail, _ = np.histogram(ec_phi[~ok & deep], h_bins)
    h_promo_ok, _ = np.histogram(ec_phi[ok & promo], h_bins)
    h_deep_ok, _ = np.histogram(ec_phi[ok & deep], h_bins)
    h_deep = h_deep_fail + h_deep_ok
    h_promo = h_promo_fail + h_promo_ok
    
    h_bins1 = 0.5 * (h_bins[:-1] + h_bins[1:])
    
    # Create dual-axis plot
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twinx()
    
    # Left axis: stacked bar chart
    width = h_bins[1] - h_bins[0]
    ax1.bar(h_bins1, h_deep_fail, width=width, label='Deep-dive Targets')
    ax1.bar(h_bins1, h_promo_fail, width=width, bottom=h_deep_fail, label='Promoted Targets')
    
    # Right axis: failure probability lines
    # Avoid division by zero
    deep_fail_rate = np.divide(h_deep_fail, h_deep, out=np.zeros_like(h_deep_fail, dtype=float), 
                               where=h_deep != 0)
    promo_fail_rate = np.divide(h_promo_fail, h_promo, out=np.zeros_like(h_promo_fail, dtype=float),
                                where=h_promo != 0)
    ax2.plot(h_bins1, deep_fail_rate, linewidth=2)
    ax2.plot(h_bins1, promo_fail_rate, linewidth=2)
    
    # Style the plot
    ax1.grid(True)
    title1 = 'Failed Earth Characterizations vs. Lambertian Phi (Whole Ensemble)'
    ax1.set_title(f'{title_x}\n{title1}', fontsize=11*1.1, fontweight='bold')
    ax1.set_xlabel('Lambertian Reflectance Phi', **tprops)
    ax1.set_ylabel('Number of Failed Characterizations [count]', **tprops)
    ax2.set_ylabel('Failure Probability', fontsize=13, **tprops)
    ax1.tick_params(labelsize=13)
    ax2.tick_params(labelsize=13)
    ax1.legend()
    
    write_plots(fig, 'earth-char-hist-phi')
    plt.close(fig)
    
    ################################################################
    # Histogram of Log Phi
    ################################################################
    
    fig, ax = plt.subplots(figsize=(8.5, 5))
    
    # Brighten jet colormap
    plt.set_cmap('jet')
    
    delta = 0.4
    h_bins_log = np.arange(-6, 0 + delta, delta)
    h_promo, _ = np.histogram(2.5 * np.log10(ec_phi[~ok & promo]), 
                              bins=np.concatenate([[-np.inf], h_bins_log]))
    h_deep, _ = np.histogram(2.5 * np.log10(ec_phi[~ok & deep]),
                            bins=np.concatenate([[-np.inf], h_bins_log]))
    
    h_bins1 = 0.5 * (h_bins_log[:-1] + h_bins_log[1:])
    h_bins1 = np.concatenate([[h_bins1[0] - delta], h_bins1])  # add leftmost bin
    
    width = delta
    ax.bar(h_bins1, h_deep, width=width, label='Deep-dive Targets')
    ax.bar(h_bins1, h_promo, width=width, bottom=h_deep, label='Promoted Targets')
    
    # Style the plot
    ax.grid(True)
    title1 = 'Failed Earth Characterizations vs. Magnitude Difference (Whole Ensemble)'
    ax.set_title(f'{title_x}\n{title1}', fontsize=11*1.1, fontweight='bold')
    ax.set_xlabel('Lambertian Reflectance Phi, as Magnitude (2.5 log$_{{10}}$ ?)', **tprops)
    ax.set_ylabel('Number of Failed Characterizations [count]', **tprops)
    ax.legend()
    ax.tick_params(labelsize=13)
    
    write_plots(fig, 'earth-char-hist-phi-log')
    plt.close(fig)

    return tracker.get_files()


def main():
    """
    Command-line interface for plot_drm_earth_chars
    """
    parser = argparse.ArgumentParser(
        description='Plot attempted earth characterizations in a DRM-set',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
    python plot_drm_earth_chars.py "data/%s.%s" "output/det-%s.%s"
    
The first argument is the source template with two %%s placeholders that will be
filled with ("info", "csv") and ("chars", "csv") to locate input files.

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

    # Load CSV data and run the plotting function
    plot_data = cs.load_csv_files(args.src_tmpl, ['earth-char-list'])
    try:
        plot_drm_earth_chars(reduce_info, plot_data, args.dest_tmpl, mode)
    except Exception as e:
        print(f"{PROGNAME}: Fatal: Unexpected error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
