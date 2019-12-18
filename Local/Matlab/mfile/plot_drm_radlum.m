function [h_bb,h_eb]=plot_drm_radlum(dest_tmpl, mode, t_info, t_radlum, t_earth)
%plot_drm_radlum	plot radius/luminosity summary for a drm-set
% 
% [h_bb,h_eb]=plot_drm_radlum(dest_tmpl, mode, t_info, t_radlum, t_earth)
% * Radius/luminosoty plots of detections.
% * If dest_tmpl is empty, no files are written.
% 
% Inputs:
%   string dest_tmpl -- filename template for outputs
%   string mode -- empty is OK
%   table t_info
%   table t_radlum -- binned detection-count table
%   table t_earth -- exo-Earth detection-count table
% 
% Outputs:
%   graphics_handle h_bb -- bar plot handles
%   graphics_handle h_eb -- error-bar plot handles
%   (files to disk)
% 
% See Also:  

% 
% Error checking
% 
if all(nargin  ~= [5]), error ('Bad input arg number'); end
% if all(nargout ~= [0 1]), error ('Bad output arg number'); end  
% allow passing an empty metadata table as []
if ismatrix(t_info) && isempty(t_info),
    t_info = cell2table(cell(0));
end;

% skip, unless mode.op contains our name or a *
if length(strfind(mode.op, '*')) == 0 && ...
   length(strfind(mode.op, 'radlum')) == 0,
    fprintf('Radius/luminosity plots: skipping, as directed.\n');
    return;
end
   
%%
%% Common data and functions
%% 

% number of radii obtained by the bin boundaries
N_rad = length(unique(t_radlum{:, 'h_Rp_lo'}));
% number of luminosities
%   obtain from the number of total bins, divided by N_rad
%   cannot use #uniques because L bin boundaries vary by radius
N_lum = length(t_radlum{:, 'h_L_lo'}) / N_rad;

% does not seem useful any more
% fprintf('%s: Got %d radius bins X %d luminosity bins.\n', mfilename, N_rad, N_lum);

% x-axis with larger inter-radius spacing 
% e.g.: [1 2 3 4.2 5.2 6.2 7.4 8.4 9.4]
x_bin_xtra = 0.29;
x_bin = reshape(bsxfun(@plus, reshape(1:N_lum*N_rad, N_lum, N_rad), ...
                       [0:N_rad-1]*x_bin_xtra), 1, []);
x_bin_plus = [-x_bin_xtra x_bin]; % includes earths

% bar color (RGB triple) for each luminosity: hot, warm, cold
%  (Stark plot uses Colors{2} =~ [0.2 0.2 0.6])
% Colors = {[0.7 0.2 0], [0.9 0.6 0.2], [0.5 0.6 1.0]};
Colors = {[0.8667 0.1843 .1020], [0.3686 0.6706 0.9608], [0.8275 0.9412 1.0]};
% for axis labels
Lumens = {'Hot', 'Warm', 'Cold'};
Lumens_plus = cat(2, {'Earth'}, repmat(Lumens, [1 N_rad]));

% error-bar color
color_eb = [0.3 0.3 0.3];
% error-bar style
style_eb = {'LineWidth', 2.5, 'CapSize', 12, 'Marker', 'none', 'Color', color_eb};

% text-block style
style_block = {'Units', 'normalized', 'FontSize', 12, 'EdgeColor', 'k', 'BackgroundColor', [0.8 0.8 0.8]};

% Helper matrix for stacked bar plots
% collapses a set of counts down
%   makes a 9x9 diagonal into a 9x3 set of stacked plots
Collapse = repmat(eye(N_lum), [N_rad 1]); 


% Inner function: rectangle underlay
% this function plots a series of vertical stripes in radius/luminosity
% plots so that planets of the same radius class are visually grouped.
    function ax_rect=place_rect_underlay(ax, N_rad, N_lum, x_bin, x_bin_xtra)

    % create an axis for the stripe underlay
    % it precisely overlaps the main axis, with the same plot units
    ax_rect = axes('Units', 'normalized', 'Position', get(ax, 'Position'));
    % set axis units
    set(ax_rect, 'XLim', get(ax, 'XLim'));
    set(ax_rect, 'YLim', get(ax, 'YLim'));
    % turn off axis adornments
    set(ax_rect, 'Color', 'none'); % make the axis itself transparent
    set(ax_rect, 'Visible', 'off', 'XTick', [], 'YTick', []);
    % put this axis on bottom
    uistack(ax_rect, 'bottom'); 

    % save plot height and width, for the rectangles below
    xlim = get(ax, 'XLim');
    ylim = get(ax, 'YLim');
    % add Earth (0) plus N_rad (1..N_rad) rectangles 
    % - one for each radius class
    for k1 = 0:N_rad,
        % rectangle is between x1 and x2, and top-to-bottom
        if k1 == 0,
            % special-case Earth
            x1 = x_bin_plus(1) - 0.5 - 0.3; % 0.5 = half-bin; 0.3 =
                                            % empirical bar width
            x2 = x_bin(1) - 0.5*(x_bin_xtra + 1); % end: left of first x_bin
            fc = [0.7569 0.8863 0.7765]; % pale green
        else,
            fc = [0.85 0.85 0.85]; % pale gray
            % index into bins: 1, N_lum+1, ...
            k2 = (k1 - 1) * N_lum + 1; 
            % we pad by 1 + x_bin_xtra between radius groups
            x1 = x_bin(k2)  - 0.5*(x_bin_xtra + 1);
            x2 = x1 + N_lum + x_bin_xtra;
            %x2 = x1 + N_lum + 0*0.5*(x_bin_xtra + 1); - 0.5*(1-x_bin_xtra)
            x2 = min(x2, xlim(2)); % cap rightmost edge at plot end
        end;
        % plot the rectangle itself
        x_gap = 0.11; % (half the) gap between adjacent rectangular panels
        % the rectangle is a child of the new axis
        % rectangle position is given as four numbers in two pairs: 
        %    [lower-left corner, width, height]
        hr = rectangle(ax_rect, ...
                       'Position',  [x1+x_gap ylim(1) (x2-x1)-2*x_gap ylim(2)-ylim(1)], ... 
                       'FaceColor', fc, ... 
                       'LineStyle', 'none');
    end;
    end


% Inner function:
%   write the current figure to various files
%   uses the dest_tmpl mfile variable
    function write_plots(dest_name)

    % write the plot out
    if ~isempty(dest_tmpl),
        % all_exts = {'png', 'pdf'};
        % pdf's don't render correctly with planet overlays
        all_exts = {'png'}; 
        for ext = all_exts,
            fn_gfx = sprintf(dest_tmpl, dest_name, ext{1});
            fprintf(1, '\tExporting to %s\n', fn_gfx);
            export_fig('-transparent', '-m2', fn_gfx);
        end;
    end
    end

% Inner function:
%   Set up plot/axis styles, title, axis labels.
%   Uses a few global mfile variables, such as t_info and x_bin.
%   Alters ax_root (typically gca)
%   Requires obs_type, a text label
    function style_bar_plot(ax_root, obs_type)
    
    ax_root = gca; % the main axis
    set(ax_root, 'YGrid', 'on');
    % XTicks - labels will repeat for as many ticks are present
    set(ax_root, 'XTick', x_bin_plus);
    set(ax_root, 'XTickLabel', Lumens_plus);
    set(ax_root, 'XLim', [x_bin_plus(1)-1 x_bin(end)+1]); % -xbin for exo-Earths
    set(ax_root, 'YLim', max(0, get(ax_root, 'YLim'))); % errorbars can go <0

    title1 = sprintf('%s vs. Insolation and Planet Radius', obs_type);
    title2 = plot_make_title(t_info);
    % title: prevent special interpretation of _
    title({title2, title1}, 'Interpreter', 'none');
    xlabel('Planet Type and Size', 'FontWeight', 'bold');
    ylabel(sprintf('%s', obs_type), 'FontWeight', 'bold');
    set(ax_root,'TitleFontSizeMultiplier', 1.1);
    set(ax_root, 'FontSize', 13);
    end % inner function

%%
%% Overall Planet Population
%% 

% data, for convenience
H_counts     = t_radlum{:, 'h_RpL_population_mean'};
H_stds       = t_radlum{:, 'h_RpL_population_std'};

% Basic plot method: use a bar-plot with overlaid error bars
% separates H_counts into N_rad groups which will have repeating colors
H_grouped = diag(H_counts) * Collapse;

% best to set up the size (aspect ratio is last two numbers)
clf;
set(gcf, 'Position', [100 100 950 500]);
   
% basic bar plot
h_bb = bar(x_bin, H_grouped, 'Stacked');
% color the bars according to luminosity
for k = 1:length(h_bb),
    set(h_bb(k), 'FaceColor', Colors{1+rem(k-1, length(Colors))});
end;

hold on

% exo-earth barplot and errorbar
%   x range consists of only a scalar
ct_e = t_earth{:,'exoE_population_mean'};
eb_e = t_earth{:,'exoE_population_std'};
bar([-x_bin_xtra], [ct_e], 'FaceColor', [0.1 0.7 0.3]);
h_ebe = errorbar(-x_bin_xtra, ct_e, eb_e, '.');
set(h_ebe, style_eb{:});

% first-stack errorbars
h_eb = errorbar(x_bin, H_counts, H_stds, '.');
set(h_eb, style_eb{:});

% Establish plot/axis styles + labels
style_bar_plot(gca, 'Planet Occurrence Rate');

% explanatory legend
block = {'Counting Over All Planets in Target List', ...
         'Normalization: Planets of that Type, per Star', ...
         'Earth Bar is Observed \eta Earth', 'Error Bar: \pm 1 sigma'};
text(0.01, 0.92, block, style_block{:});

% planet + rectangle underlays
ax_main = gca;
ax_ulay = plot_drm_planet_overlay(ax_main, mode);
ax_rlay = place_rect_underlay(ax_main, N_rad, N_lum, x_bin, x_bin_xtra);

% write the plots
write_plots('radlum-population');


%%
%% Detections - unique
%% 

% data, for convenience
H_counts     = t_radlum{:, 'h_RpL_det_main_mean'};
H_stds       = t_radlum{:, 'h_RpL_det_main_std'};
H_counts_alt = t_radlum{:, 'h_RpL_det_alt_mean'};
H_stds_alt   = t_radlum{:, 'h_RpL_det_alt_std'};
% pooled standard error of the sum: std = sqrt(std^2 + std^2)
H_stds_alt2  = hypot(H_stds, H_stds_alt);

% Basic plot method: use a bar-plot with overlaid error bars
% separates H_counts into N_rad groups which will have repeating colors
H_grouped = diag(H_counts) * Collapse;
H_grouped_alt = diag(H_counts_alt) * Collapse;

% best to set up the size (aspect ratio is last two numbers)
clf;
set(gcf, 'Position', [100 100 950 500]);
   
% basic bar plot
h_bb = bar(x_bin, [H_grouped, H_grouped_alt], 'Stacked');
% color the bars according to luminosity
for k = 1:length(h_bb),
    set(h_bb(k), 'FaceColor', Colors{1+rem(k-1, length(Colors))});
end;

hold on

% exo-earth barplot and errorbar
%   the "nan" construct is a work-around to allow stacked bars when the
%   x range consists of only a scalar
ct_e = [t_earth{:,'exoE_det_main_mean'}, t_earth{:,'exoE_det_alt_mean'}];
eb_e = hypot(t_earth{:,'exoE_det_main_std'}, t_earth{:,'exoE_det_alt_std'});
bar([nan -x_bin_xtra], [nan(1,2); ct_e], ...
    'Stacked', 'FaceColor', [0.1 0.7 0.3]);
h_ebe = errorbar(-x_bin_xtra, sum(ct_e), eb_e, '.');
set(h_ebe, style_eb{:});

% first-stack errorbars
h_eb = errorbar(x_bin, H_counts, H_stds, '.');
set(h_eb, style_eb{:});
% second-stack errorbars (should be optional)
h_eb_alt = errorbar(x_bin-0.2, H_counts+H_counts_alt, H_stds_alt, '.');
set(h_eb_alt, style_eb{:});
% pooled errorbars (again, should be optional)
h_eb_alt2 = errorbar(x_bin+0.2, H_counts+H_counts_alt, H_stds_alt2, '.');
set(h_eb_alt2, style_eb{:});

% Establish plot/axis styles + labels
style_bar_plot(gca, 'Unique Detections');

% explanatory legend
block = {'Lower Segment: Blue Detection Mode', 'Upper Segment: Combo Detection Mode', 'Error Bar: \pm 1 sigma'};
text(0.01, 0.92, block, style_block{:});

% planet + rectangle underlays
ax_main = gca;
ax_ulay = plot_drm_planet_overlay(ax_main, mode);
ax_rlay = place_rect_underlay(ax_main, N_rad, N_lum, x_bin, x_bin_xtra);

% write the plots
write_plots('radlum-det');

%%
%% Detections - all
%% 

% (copy-paste from above, alas)

% data, for convenience
H_counts     = t_radlum{:, 'h_RpL_xdet_main_mean'};
H_stds       = t_radlum{:, 'h_RpL_xdet_main_std'};
H_counts_alt = t_radlum{:, 'h_RpL_xdet_alt_mean'};
H_stds_alt   = t_radlum{:, 'h_RpL_xdet_alt_std'};
% pooled standard error of the sum: std = sqrt(std^2 + std^2)
H_stds_alt2  = hypot(H_stds, H_stds_alt);

% Basic plot method: use a bar-plot with overlaid error bars
% separates H_counts into N_rad groups which will have repeating colors
H_grouped = diag(H_counts) * Collapse;
H_grouped_alt = diag(H_counts_alt) * Collapse;

% best to set up the size (aspect ratio is last two numbers)
clf;
set(gcf, 'Position', [100 100 950 500]);
   
% basic bar plot
h_bb = bar(x_bin, [H_grouped, H_grouped_alt], 'Stacked');
% color the bars according to luminosity
for k = 1:length(h_bb),
    set(h_bb(k), 'FaceColor', Colors{1+rem(k-1, length(Colors))});
end;

hold on

% exo-earth barplot and errorbar
%   the "nan" construct is a work-around to allow stacked bars when the
%   x range consists of only a scalar
ct_e = [t_earth{:,'exoE_xdet_main_mean'}, t_earth{:,'exoE_xdet_alt_mean'}];
eb_e = hypot(t_earth{:,'exoE_xdet_main_std'}, t_earth{:,'exoE_xdet_alt_std'});
bar([nan -x_bin_xtra], [nan(1,2); ct_e], ...
    'Stacked', 'FaceColor', [0.1 0.7 0.3]);
h_ebe = errorbar(-x_bin_xtra, sum(ct_e), eb_e, '.');
set(h_ebe, style_eb{:});

% first-stack errorbars
h_eb = errorbar(x_bin, H_counts, H_stds, '.');
set(h_eb, style_eb{:});
% second-stack errorbars (should be optional)
h_eb_alt = errorbar(x_bin-0.2, H_counts+H_counts_alt, H_stds_alt, '.');
set(h_eb_alt, style_eb{:});
% pooled errorbars (again, should be optional)
h_eb_alt2 = errorbar(x_bin+0.2, H_counts+H_counts_alt, H_stds_alt2, '.');
set(h_eb_alt2, style_eb{:});

% Establish plot/axis styles + labels
style_bar_plot(gca, 'All Detections');

% explanatory legend
block = {'Lower Segment: Blue Detection Mode', 'Upper Segment: Combo Detection Mode', 'Error Bar: \pm 1 sigma'};
text(0.01, 0.92, block, style_block{:});

% planet + rectangle underlays
ax_main = gca;
ax_ulay = plot_drm_planet_overlay(ax_main, mode);
ax_rlay = place_rect_underlay(ax_main, N_rad, N_lum, x_bin, x_bin_xtra);

% write the plots
write_plots('radlum-det-all');

%%
%% Characterizations
%% 

%% Controls the characterization plots we make.
% Variable format in each row below:
%  (1) segment of variable name in the t_radlum variable
%  (2) word for title of plot
%  (3) full/partial w/ stacked bars (dual) vs. solo (no stacked bars)
%  (4) output filename component
% We chose the "union" version as the default version, so we leave its
% descriptive title-word blank, to signal that this qualifier is omitted 
% in the plot title as shown.
char_plot_menu = {...
    {'_',        'Legacy', 'dual', 'radlum-char'}, ...
    {'strict_',  'Strict', 'solo', 'radlum-char-strict'}, ...
    {'_union_',  '',       'dual', 'radlum-char-union'}, ...
    {'_blue_',   'Blue',   'dual', 'radlum-char-blue'}, ...
    {'_red_',    'Red',    'dual', 'radlum-char-red'}, ...
    };

for plot_num = 1:length(char_plot_menu),
    if length(char_plot_menu{plot_num}) ~= 4, 
        error(sprintf('Plot selector %d has too few elements', plot_num));
    end;
    plot_prop      = char_plot_menu{plot_num}{1};
    plot_title_raw = char_plot_menu{plot_num}{2};
    plot_mode      = char_plot_menu{plot_num}{3};
    plot_file      = char_plot_menu{plot_num}{4};

    if isempty(plot_title_raw),
        plot_title = '';
    else,
        plot_title = sprintf(' [%s]', plot_title_raw);
    end;

    % not all bands will necessarily be present in the t_radlum{} table
    if strcmp(plot_mode, 'dual'),
        plot_prop_root = '_full';
    else,
        plot_prop_root = '_'; % (strict)
    end
    test_name = sprintf('h_RpL_char%s%smean', plot_prop_root, plot_prop);
    if ~any(strcmp(test_name, t_radlum.Properties.VariableNames)),
        continue;
    end

    % data, for convenience
    %   full and partial characterizations
    %   it is just coincidence that there are again two types of stacked bars
    H_counts     = t_radlum{:, sprintf('h_RpL_char%s%smean', plot_prop_root, plot_prop)};
    H_stds       = t_radlum{:, sprintf('h_RpL_char%s%sstd',  plot_prop_root, plot_prop)};
    if strcmp(plot_mode, 'dual'),
        H_counts_alt = t_radlum{:, sprintf('h_RpL_char_part%smean', plot_prop)};
        H_stds_alt   = t_radlum{:, sprintf('h_RpL_char_part%sstd',  plot_prop)};
        % pooled standard error of the sum: std = sqrt(std^2 + std^2)
        H_stds_alt2  = hypot(H_stds, H_stds_alt);
    else,
        H_counts_alt = [];
        H_stds_alt = [];
        H_stds_alt2 = [];
    end;

    % Basic plot method: use a bar-plot with overlaid error bars
    % separates H_counts into N_rad groups which will have repeating colors
    H_grouped = diag(H_counts) * Collapse;
    if ~isempty(H_counts_alt),
        H_grouped_alt = diag(H_counts_alt) * Collapse;
    else,
        H_grouped_alt = [];
    end;

    % best to set up the size (aspect ratio is last two numbers)
    clf; set(gcf, 'Position', [100 100 950 500]);
    
    % basic bar plot -- stack of one bar in solo mode
    h_bb = bar(x_bin, [H_grouped, H_grouped_alt], 'Stacked');
    % color the bars according to luminosity
    for k = 1:length(h_bb),
        set(h_bb(k), 'FaceColor', Colors{1+rem(k-1, length(Colors))});
    end;

    hold on

    % exo-earth barplot and errorbar
    %   the "nan" construct is a work-around to allow stacked bars when the
    %   x range consists of only a scalar
    if strcmp(plot_mode, 'dual'),
        ct_e = [t_earth{:, sprintf('exoE_char_full%smean', plot_prop)}, ...
                t_earth{:, sprintf('exoE_char_part%smean', plot_prop)}];
        eb_e = hypot(t_earth{:, sprintf('exoE_char_full%sstd', plot_prop)},...
                     t_earth{:, sprintf('exoE_char_part%sstd', plot_prop)});
        bar([nan -x_bin_xtra], [nan(1,2); ct_e], ...
            'Stacked', 'FaceColor', [0.1 0.7 0.3]);
        h_ebe = errorbar(-x_bin_xtra, sum(ct_e), eb_e, '.');
    else,
        ct_e = t_earth{:, sprintf('exoE_char_%smean', plot_prop)};
        eb_e = t_earth{:, sprintf('exoE_char_%sstd',  plot_prop)};
        bar(-x_bin_xtra, ct_e, ...
            'FaceColor', [0.1 0.7 0.3]);
        h_ebe = errorbar(-x_bin_xtra, ct_e, eb_e, '.');
    end;
    set(h_ebe, style_eb{:});

    % first-stack errorbars
    h_eb = errorbar(x_bin, H_counts, H_stds, '.');
    set(h_eb, style_eb{:});
    if ~isempty(H_counts_alt),
        % second-stack errorbars (optional)
        h_eb_alt = errorbar(x_bin-0.2, H_counts+H_counts_alt, H_stds_alt, '.');
        set(h_eb_alt, style_eb{:});
        % combo errorbars (optional)
        h_eb_alt2 = errorbar(x_bin+0.2, H_counts+H_counts_alt, H_stds_alt2, '.');
        set(h_eb_alt2, style_eb{:});
    end;

    % Establish plot/axis styles + labels
    % This will apply to the left axis
    style_bar_plot(gca, sprintf('Unique Characterizations%s', plot_title));

    % explanatory legend
    block_top = sprintf('Unique Characterizations%s', plot_title);
    block_eb = '  Error Bar: \pm 1 sigma';
    if ~isempty(H_counts_alt),
        block = {block_top, ...
                 '  Lower Segment: Full Characterization', ...
                 '  Upper Segment: Partial Characterization', ...
                 block_eb};
    else,
        block = {block_top, ...
                 block_eb};
    end;
    text(0.01, 0.90, block, style_block{:});

    % planet + rectangle underlays
    ax_main = gca;
    ax_ulay = plot_drm_planet_overlay(ax_main, mode);
    ax_rlay = place_rect_underlay(ax_main, N_rad, N_lum, x_bin, x_bin_xtra);

    % write the plots
    write_plots(plot_file);

end;

%%
%% Characterization SNR
%% 
  
%% Controls the characterization plots we make.
% Variable format in each row below:
%  (1) segment of variable name in the t_radlum variable
%  (2) word for title of plot
%  (3) output filename component
% We chose the "union" version as the default version, so we leave its
% descriptive title-word blank, to signal that this qualifier is omitted 
% in the plot title as shown.
snr_plot_menu = {...
    {'_',       'Legacy', 'radlum-char-snr'}, ...
    {'_union_', '',       'radlum-char-snr-union'}, ...
    {'_blue_',  'Blue',   'radlum-char-snr-blue'}, ...
    {'_red_',   'Red',    'radlum-char-snr-red'}, ...
    };

for plot_num = 1:length(snr_plot_menu),
    if length(snr_plot_menu{plot_num}) ~= 3, 
        error(sprintf('Plot selector %d has too few elements', plot_num));
    end;
    plot_prop      = snr_plot_menu{plot_num}{1};
    plot_title_raw = snr_plot_menu{plot_num}{2};
    plot_file      = snr_plot_menu{plot_num}{3};

    if isempty(plot_title_raw),
        plot_title = '';
    else,
        plot_title = sprintf(' [%s]', plot_title_raw);
    end;

    % not all bands will necessarily be present in the table
    test_name = sprintf('h_RpL_char_full%smean', plot_prop);
    if ~any(strcmp(test_name, t_radlum.Properties.VariableNames)),
        continue;
    end

    clf; set(gcf, 'Position', [100 100 950 500]);
   

    % Characterization SNR plot
    %   note, different marker style 
    % quantile error bar
    h_eb_right = errorbar(x_bin, ...
                          t_radlum{:,sprintf('h_RpL_char_snr%sq50', plot_prop)}, ...
                          t_radlum{:,sprintf('h_RpL_char_snr%sq50', plot_prop)} - t_radlum{:,sprintf('h_RpL_char_snr%sq25', plot_prop)}, ...
                          t_radlum{:,sprintf('h_RpL_char_snr%sq75', plot_prop)} - t_radlum{:,sprintf('h_RpL_char_snr%sq50', plot_prop)}, ...
                          'o', 'MarkerSize', 10);
    % add Earth
    hold on;
    h_eb_rightE = errorbar(-x_bin_xtra, ...
                           t_earth{:,sprintf('exoE_char_snr_q50', plot_prop)}, ...
                           t_earth{:,sprintf('exoE_char_snr_q50', plot_prop)} - t_earth{:,sprintf('exoE_char_snr_q25', plot_prop)}, ...
                           t_earth{:,sprintf('exoE_char_snr_q75', plot_prop)} - t_earth{:,sprintf('exoE_char_snr_q50', plot_prop)}, ...
                           'o', 'MarkerSize', 10);

    % style, labels
    style_eb2 = {'LineWidth', 2.5, 'CapSize', 12, 'Color', color_eb};
    set(h_eb_right,  style_eb2{:});
    set(h_eb_rightE, style_eb2{:});
    % Limit plot y-range
    % (A) SNR cannot be negative, but error bars sometimes go < 0
    % (B) SNR upper YLim can sometimes be much too high, because SNR std is very high
    top_snr = 1.1 * max(t_radlum{:,sprintf('h_RpL_char_snr%smean', plot_prop)});
    ylim_snr = get(gca, 'YLim');
    set(gca, 'YLim', [max(0, ylim_snr(1)) min(ylim_snr(2), top_snr)]);

    % Establish plot/axis styles + labels
    % This will apply to the left axis
    style_bar_plot(gca, sprintf('Characterization SNR%s', plot_title));

    % explanatory legend
    block = {sprintf('Characterization SNR%s', plot_title), ...
             '  Circles: Median', ...
             '  Error Bars: Q25, Q75'};
    text(0.01, 0.90, block, style_block{:});

    % planet + rectangle underlays
    ax_main = gca;
    ax_ulay = plot_drm_planet_overlay(ax_main, mode);
    ax_rlay = place_rect_underlay(ax_main, N_rad, N_lum, x_bin, x_bin_xtra);

    % put the axis grid on top of the above annotations
    ax_main.Layer = 'top';

    % write the plots
    write_plots(plot_file);

end

return

end % main function
