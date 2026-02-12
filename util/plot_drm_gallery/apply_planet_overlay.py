#!/usr/bin/env python

# Still here for reference
# This is the original converted file - this was inserted into
# plot_drm_radlum.py and modified there

"""
Plot a three-planet decorative underlay

Plot a decorative underlay onto ax_dest (typically current axes),
returning the axis object of the underlay.
The mode is currently unused.
"""

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from pathlib import Path
import numpy as np
import os


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
        ctrline = 0.72  # centerline of the left-to-right line of planets
        ppos = [
            [0.24, ctrline - 0.033/2, 0.0333, 0.0333],
            [0.35, ctrline - 0.1167/2, 0.1167, 0.1167],
            [0.42, ctrline - 0.233/2, 0.2333, 0.2222],
            [0.48, ctrline - 0.40/2, 0.40, 0.40],
            [0.63, ctrline - 0.73*j_pre/2 + 0.01, 0.73*j_pre, 0.73*j_pre]
        ]
    
    # Make planet underlay
    ax_list = []
    
    # Kill white background on current axes, so underlay is visible
    ax_dest.patch.set_facecolor('none')
    
    for p, pfile_ in enumerate(pfiles):
        # FIXME: temporary hack for resources
        pfile = 'Local/Matlab/mfile' / Path(pfile_)
        # Check if file exists
        if not os.path.exists(pfile):
            print(f"Warning: Planet image file '{pfile}' not found, skipping")
            continue
        
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
            ax1 = inset_axes(ax_dest, width=f'{pos[2]*100}%', height=f'{pos[3]*100}%',
                           loc='lower left',
                           bbox_to_anchor=(pos[0], pos[1], pos[2], pos[3]),
                           bbox_transform=ax_dest.transAxes,
                           borderpad=0)
            
            # Plot the image into the axis (including transparency info)
            ax1.imshow(im, alpha=alfa, aspect='equal')
            
            # Set axis parameters
            ax1.set_facecolor('none')  # make the axis itself transparent
            ax1.set_aspect('equal')    # 1:1 aspect ratio so planet-images are circles
            
            # Turn off axis adornments
            ax1.set_xticks([])
            ax1.set_yticks([])
            ax1.axis('off')
            
            # Set zorder to put on bottom
            ax1.set_zorder(-1)
            
            # Record the underlay axis, ax1
            ax_list.append(ax1)
            
        except Exception as e:
            print(f"Warning: Could not load planet image '{pfile}': {e}")
            continue
    
    return ax_list


# Example standalone usage for testing
def main():
    """Test the planet overlay function"""
    import matplotlib.pyplot as plt
    
    # Create a simple plot to overlay planets on
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Sample data
    x = np.linspace(0, 10, 100)
    y = np.sin(x)
    ax.plot(x, y)
    
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_title('Test Plot with Planet Overlay')
    ax.grid(True)
    
    # Add planet overlay
    mode = {}
    ax_planets = plot_drm_planet_overlay(ax, mode)
    
    print(f"Created {len(ax_planets)} planet overlays")
    
    plt.tight_layout()
    plt.savefig('test_planet_overlay.png', dpi=150, bbox_inches='tight')
    plt.show()


if __name__ == '__main__':
    main()
