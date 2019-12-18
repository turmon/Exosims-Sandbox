function ax=plot_drm_planet_overlay(ax_dest, mode)
%plot_drm_planet_overlay	plot a three-planet decorative underlay
% 
% ax=plot_drm_planet_overlay(ax_dest, mode)
% * Plot a decorative underlay onto ax_dest (typically "gca"),
% returning the axis object of the underlay.
% * The mode is currently unused.
% 
% Inputs:
%   graphics_handle ax_dest -- destination axes object
%   opt struct mode = {}
% 
% Outputs:
%   graphics_handle ax -- planet axis
% 
% See Also:  

% 
% Error checking
% 
if all(nargin  ~= [1 2]), error ('Bad input arg number'); end
% if all(nargout ~= [0 1]), error ('Bad output arg number'); end  
if nargin < 2, 
    mode = struct(); 
end;

%
% Computation
% 

if false,
    %% 3-class setup (late 2017, early 2018)
    % planets: earthlike, neptune, jupiter, left-to-right
    pfiles = {'gliese-667C_c.png', 'neptune.png', 'jupiter.png'};
    % positions are x0, y0, xsize, ysize
    ppos = {[0.15 0.6 0.15 0.15], [0.35 0.6 0.3 0.3], [0.58 0.45 0.45 0.45]};
    % pixels to crop each planet (above) on each edge
    pcrops = {[0 0 0 0], [70 70 70 70], [300 80 80 400]};
elif false,
    %% 5-class setup (mid 2018)
    % planets, left-to-right: 
    % rocky sub-earth, super-earth, sub-neptune, neptune, jupiter
    pfiles = {'trappist-1_d.png', 'gliese-667C_c.png', 'gj-1214_b.png', 'neptune.png', 'jupiter.png'};
    % positions are x0, y0, xsize, ysize
    % (older, smaller planets)
    ppos = {[0.22 0.4 0.10 0.10], [0.33 0.50 0.15 0.15], [0.45 0.6 0.2 0.2], [0.53 0.62 0.3 0.3], [0.63 0.57 0.35 0.35]};
    % pixels to crop each planet on each edge (top, left, bottom, right)
    pcrops = {[25 25 25 25], [0 0 0 0], [0 0 0 0], [70 70 70 70], [300 80 80 300]};
else,
    %% 5-class setup, planets to scale (july 2018)
    % planets, left-to-right: 
    % rocky sub-earth, super-earth, sub-neptune, neptune, jupiter (NB: cropped tight)
    pfiles = {'trappist-1_d.png', 'gliese-667C_c.png', 'gj-1214_b.png', 'neptune.png', 'jupiter-tight.png'};
    % pixels to crop each planet (above) on each edge (top, left, bottom, right)
    % Note: All crops except Jupiter are to the planetary disk, so that
    % scaling is understood
    pcrops = {[25 25 25 25], [25 25 25 25], [10 10 10 10], [100 100 100 100], [300 0 0 300]};
    J_pre = 0.5; % amount the Jupiter-class image is
                 % already shrunken due to cropping - others are
                 % cropped on the planetary disk
    % Set planet positions within the enclosing axis
    % positions are x0, y0, xsize, ysize
    % size is relative to radius property of the corresponding planet class
    ctrline = 0.72; % centerline of the left-to-right line of planets
    ppos = {[0.24  ctrline-.033/2            0.0333 0.0333], ...
            [0.35  ctrline-.1167/2           0.1167 0.1167], ...
            [0.42  ctrline-0.233/2           0.2333 0.2222], ...
            [0.48  ctrline-0.40/2            0.40   0.40], ...
            [0.63  ctrline-0.73*J_pre/2+.01  0.73*J_pre 0.73*J_pre]};
end;    

%% make planet underlay
ax = [];
if true,
    % kill white background on current axes, so underlay is visible
    set(ax_dest, 'Color', 'none');
    for p = 1:length(pfiles),
        % load underlay image
        [im_orig,~,alfa_orig] = imread(pfiles{p});
        % crop image and alpha map
        pcrop = pcrops{p}; % crop info, for convenience
        im   = im_orig  (1+pcrop(1):end-pcrop(3), 1+pcrop(2):end-pcrop(4),:);
        alfa = alfa_orig(1+pcrop(1):end-pcrop(3), 1+pcrop(2):end-pcrop(4),:);
        % create an axis for the underlay
        ax1 = axes('Units', 'normalized', 'Position', ppos{p});
        % put this axis on bottom
        uistack(ax1, 'bottom'); 
        % plot the image into the axis (including transparency info)
        h_im = image(ax1, im, 'AlphaData', alfa);
        % set axis parameters
        set(ax1, 'Color', 'none'); % make the axis itself transparent
        % 1:1 aspect ratio so planets-images are circles
        axis(ax1, 'equal'); 
        % turn off axis adornments
        set(ax1, 'Visible', 'off', 'XTick', [], 'YTick', []);
        % record the underlay axis, ax1
        ax(p) = ax1;
    end;
end

return
