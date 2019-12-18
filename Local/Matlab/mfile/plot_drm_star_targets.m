function plot_drm_star_targets(dest_tmpl, mode, t_info, t_star_targ)
%plot_drm_star_targets	plot yield, etc., against star targets in a drm-set
% 
% plot_drm_star_targets(dest_tmpl, mode, t_info, t_star_targ)
% * Make, and write out, plots of various quantities pertaining to each
% star target in a set of DRMs.  
% * An example plot of this type is yield for each star, plotted against
% stellar luminosity and distance.
% 
% Inputs:
%   string dest_tmpl
%   string mode
%   table t_info
%   table t_star_targ
% 
% Outputs:
%   (to disk)
% 
% See Also:  

% 
% Error checking
% 
% if all(nargin  ~= [4]), error ('Bad input arg number'); end
% if all(nargout ~= [0 1]), error ('Bad output arg number'); end  
% allow passing an empty metadata table as []
if ismatrix(t_info) && isempty(t_info),
    t_info = cell2table(cell(0));
end;

% allow skipping these plots so we don't fail on old reductions
if ismatrix(t_star_targ) && isempty(t_star_targ),
    fprintf('Star target plots: skipping (re-run reduction?).\n');
    return;
end

% skip, unless mode.op contains our name or a *
if length(strfind(mode.op, '*')) == 0 && ...
   length(strfind(mode.op, 'perstar')) == 0,
    fprintf('Star target plots: skipping, as directed.\n');
    return;
end

% 
% Common data
% 

% color map
%cmap = flipud(parula(256));
cmap0 = jet(280); cmap = cmap0(end-256:end-1,:); % initial blue too dark
% 'scatter' dot size
dot_size = 400;

% clip target plots at 30 pc
dist_axis_limit = [0.0 30.0];

% experiment name, or empty
%if any(strcmp('experiment', fieldnames(t_info))),
%    ExpName = t_info{1,'experiment'}{1};
%else,
%    ExpName = '';
%end

% ensemble size
%EnsSize = t_info{1,'ensemble_size'};

% file extensions to write
%ExtList = {'png', 'pdf'};
ExtList = {'png'};

% color for un-observed stars
unseen_color = [1 1 1]*0.7;

% 
% Common functions
% 

%% Inner function:
%   Set up plot/axis styles, title, axis labels.
%   Uses a few global mfile variables, such as t_info and x_bin.
%   Alters ax_root (typically gca)
%   Requires obs_name, a text label, and obs_unit, its unit
    function style_scatter_plot(ax_root, obs_name, obs_unit)

    %  Plot/axis styles
    title1 = sprintf('Luminosity vs. Distance, Shaded by %s', obs_name);
    title2 = plot_make_title(t_info);

    % title: prevent special interpretation of _
    title({title2, title1}, 'Interpreter', 'none');
    xlabel('Distance [pc]', 'FontWeight', 'bold')
    ylabel('Luminosity [L_{sun}]', 'FontWeight', 'bold')
    set(ax_root, 'TitleFontSizeMultiplier', 1.1)
    set(ax_root, 'FontSize', 13);
    set(ax_root, 'YScale', 'log'); % semilog
    set(ax_root, 'YLimSpec', 'Tight'); % tight on Y
    xlim(dist_axis_limit); % clip distance
    grid on;
    colormap(cmap);
    h = colorbar;
    ylabel_string = sprintf('%s [%s]', obs_name, obs_unit);
    ylabel(h, ylabel_string, 'FontWeight', 'bold');

    end % inner function
    
%% inner function: write the current plot out in whatever format(s) needed
function write_plots(file_tag);
    if ~isempty(dest_tmpl),
        for ext = ExtList,
            fn_gfx = sprintf(dest_tmpl, file_tag, ext{1});
            fprintf(1, '\tExporting to %s\n', fn_gfx);
            export_fig('-transparent', '-m2', fn_gfx);
        end;
    end
    end % inner function


%
% Computation
% 

%% Graphics setup

% set up the size (aspect ratio is last two numbers)
clf;
set(gcf, 'Position', [100 100 850 500]);
   
% zero-visit (unseen) stars are not the same as zero-yield stars
% this variable keeps track of what stars were visited
seen = (t_star_targ{:,'h_star_det_visit_mean'} > 0);
earth = (t_star_targ{:,'h_star_det_earth_cume_mean'} > 0);

%%
% The plot_menu variable controls the properties that are plotted
% in standard distance/luminosity coordinates.  Each entry is one
% plot, with a variable to plot, units, filename output, etc.
%
% All the plots look about the same, so this was the best way to 
% make the series of plots without duplicating code.
%
% Variable format in each row below:
%   'h_star_...'    -> entry in the T_star_targ table
%   passthru/@log10 -> variable transform for the color-coding
%   true/false -> put on the "we visited this" overlay (for Earths)
%   'Title'    -> plot title
%   'unit'     -> plot color-coding label (x/y axis are always dist/lum)
%   'perstar_...' -> output filename component
%
% Plots we make:
%   * Yields (2x2x2 = 8):
%     {Detection, Char.} X {All Planets, Earthlike} X {total, unique}
%   * Integration Time (2): Detection, Characterization
%   * Yield-per-time = "Value", titled as "Rank" (2x2 = 4):
%     {Detection, Char.} X {All Planets, Earthlike}

% dummy function to pass matrix through unchanged
passthru = @(x) x;

plot_menu = {...
    {'h_star_det_plan_cume_mean',  passthru, false, ...
         'Mean Total Detections',  'count', ...
         'perstar-det-allplan-cume'}, ...
    {'h_star_det_plan_uniq_mean',  passthru, false, ...
         'Mean Unique Detections',  'count', ...
         'perstar-det-allplan-uniq'}, ...
    {'h_star_det_earth_cume_mean', passthru, true, ...
         'Mean Total Detections: Earth', 'count', ...
         'perstar-det-earth-cume'}, ...
    {'h_star_det_earth_uniq_mean', passthru, true, ...
         'Mean Unique Detections: Earth', 'count', ...
         'perstar-det-earth-uniq'}, ...
    {'h_star_char_plan_cume_mean', passthru, false, ...
         'Mean Total Characterizations', 'count', ...
         'perstar-char-allplan-cume'}, ...
    {'h_star_char_plan_uniq_mean', passthru, false, ...
         'Mean Unique Characterizations', 'count', ...
         'perstar-char-allplan-uniq'}, ...
    {'h_star_char_earth_cume_mean', passthru, true, ...
         'Mean Total Characterizations: Earth', 'count', ...
         'perstar-char-earth-cume'}, ...
    {'h_star_char_earth_uniq_mean', passthru, true, ...
         'Mean Unique Characterizations: Earth', 'count', ...
         'perstar-char-earth-uniq'}, ...
    {'h_star_det_tInt_mean',  passthru, false, ...
         'Mean Integration Time (Det.)',  'day', ...
         'perstar-det-t-int'}, ...
    {'h_star_char_tInt_mean', passthru, false, ...
         'Mean Integration Time (Char.)', 'day', ...
         'perstar-char-t-int'}, ...
    {'h_star_det_plan_value_mean', @log10, false, ...
         'Detection Rank', 'log_{10} count/day', ...
         'perstar-det-allplan-rank'}, ...
    {'h_star_det_plan_frac_mean', passthru, false, ...
         'Planets Detected/Planets Present', 'count/count', ...
         'perstar-det-allplan-frac'}, ...
    {'h_star_char_plan_value_mean', @log10, false, ...
         'Characterization Rank', 'log_{10} count/day', ...
         'perstar-char-allplan-rank'}, ...
    {'h_star_char_plan_frac_mean', passthru, false, ...
         'Planets Characterized/Planets Present', 'count/count', ...
         'perstar-char-allplan-frac'}, ...
    {'h_star_det_earth_value_mean', @log10, true, ...
         'Earth Detection Rank', 'log_{10} count/day', ...
         'perstar-det-earth-rank'}, ...
    {'h_star_det_earth_frac_mean', passthru, true, ...
         'Earths Detected/Earths Present', 'count/count', ...
         'perstar-det-earth-frac'}, ...
    {'h_star_char_earth_value_mean', @log10, true, ...
         'Earth Characterization Rank', 'log_{10} count/day', ...
         'perstar-char-earth-rank'}, ...
    {'h_star_char_earth_frac_mean', passthru, true, ...
         'Earths Characterized/Earths Present', 'count/count', ...
         'perstar-char-earth-frac'}, ...
};


%%
% This loop runs through the table above

for plot_num = 1:length(plot_menu),
    if length(plot_menu{plot_num}) ~= 6, 
        error(sprintf('Plot selector %d has too few elements', plot_num));
    end;
    plot_prop  = plot_menu{plot_num}{1}; % property to plot
    plot_xform = plot_menu{plot_num}{2}; % transformation of the property
    plot_earth = plot_menu{plot_num}{3}; % is it a plot of Earthlike dets/chars?
    plot_title = plot_menu{plot_num}{4}; % title, text
    plot_unit  = plot_menu{plot_num}{5}; % unit of plot_prop above
    plot_file  = plot_menu{plot_num}{6}; % filename

    clf;
    % selector for the final, "on-top" scatter plot
    if plot_earth,
        selector = earth;
    else,
        selector = seen;
    end
    % un-visited in plain color, transparency, smaller dots
    scatter(t_star_targ{~seen,'star_dist'}, ...
            t_star_targ{~seen,'star_L'}, dot_size/4, ...
            unseen_color, '.', 'MarkerFaceAlpha', 0.2);
    hold on;
    % show visited stars, Earth or otherwise
    if plot_earth,
        scatter(t_star_targ{seen,'star_dist'}, ...
                t_star_targ{seen,'star_L'}, dot_size/2, ...
                [1 1 1]*0.4, '.', 'MarkerFaceAlpha', 0.2);

    end;
    % visited, in color varying by detection rate
    scatter(t_star_targ{selector,'star_dist'}, ...
            t_star_targ{selector,'star_L'}, dot_size, ...
            plot_xform(t_star_targ{selector,plot_prop}), '.');

    %  Plot/axis styles
    style_scatter_plot(gca, plot_title, plot_unit);
    write_plots(plot_file);
end;


%%
%%  Change to another plot format
%%


%% Alternative Plot: tInt vs. yield -- Detections
clf;
scatter(t_star_targ{seen, 'h_star_det_tInt_mean'}, ...
        t_star_targ{seen, 'h_star_det_plan_cume_mean'}, dot_size, ...
        t_star_targ{seen, 'h_star_det_tobs1_mean'}, '.');

set(gca, 'XScale', 'log');
%set(gca, 'YScale', 'log');

%  Plot/axis styles
title1 = 'Cumulative Integration Time vs. Detections, Shaded by First-Observation Time';
title2 = plot_make_title(t_info);

% title: prevent special interpretation of _
title({title2, title1}, 'Interpreter', 'none');
xlabel('Integration Time [d]', 'FontWeight', 'bold')
ylabel('Mean Yield [count]', 'FontWeight', 'bold')
set(gca,'TitleFontSizeMultiplier', 1.1)
set(gca, 'FontSize', 13);
grid on;
colormap(cmap);
h = colorbar;
ylabel(h, 'Time of First Obs. [day]', 'FontWeight', 'bold');

write_plots('perstar-det-tint-yield');

%% Alternative Plot: tInt vs. yield -- Characterizations
clf;
scatter(t_star_targ{seen, 'h_star_char_tInt_mean'}, ...
        t_star_targ{seen, 'h_star_char_plan_cume_mean'}, dot_size, ...
        t_star_targ{seen, 'h_star_char_tobs1_mean'}, '.');

set(gca, 'XScale', 'log');
%set(gca, 'YScale', 'log');

%  Plot/axis styles
title1 = 'Cumulative Integration Time vs. Characterizations, Shaded by First-Observation Time';
title2 = plot_make_title(t_info);

% title: prevent special interpretation of _
title({title2, title1}, 'Interpreter', 'none');
xlabel('Integration Time [d]', 'FontWeight', 'bold')
ylabel('Mean Characterizations [count]', 'FontWeight', 'bold')
set(gca,'TitleFontSizeMultiplier', 1.1)
set(gca, 'FontSize', 13);
grid on;
colormap(cmap);
h = colorbar;
ylabel(h, 'Time of First Obs. [day]', 'FontWeight', 'bold');

write_plots('perstar-char-tint-yield');


return
end
