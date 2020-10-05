function plot_drm_event_counts(dest_tmpl, mode, t_info, t_counts, t_earth_counts)
%plot_drm_event_counts	plot event counts in a drm-set
% 
% plot_drm_event_counts(dest_tmpl, mode, t_info, t_counts, t_earth_counts)
% * Plots of event-counts such as number of slews, characterizations, 
% detections, averaged over an ensemble.
% 
% Inputs:
%   string dest_tmpl
%   string mode
%   table t_info
%   table t_counts
%   table t_earth_counts
% 
% Outputs:
%   (to disk)
% 
% See Also:  

% 
% Error checking
% 
% if all(nargin  ~= [3]), error ('Bad input arg number'); end
% if all(nargout ~= [0 1]), error ('Bad output arg number'); end  
% allow passing an empty metadata table as []
if ismatrix(t_info) && isempty(t_info),
    t_info = cell2table(cell(0));
end;

% skip, unless mode.op contains our name or a *
if length(strfind(mode.op, '*')) == 0 && ...
   length(strfind(mode.op, 'event-counts')) == 0,
    fprintf('Event-count plots: skipping, as directed.\n');
    return;
end

% file extensions to write
%ExtList = {'png', 'pdf'};
ExtList = {'png'};

%
% Utilities
% 

% Inner function:
%   Set up plot/axis styles, title, axis labels.
%   Uses global mfile variables, such as t_info
%   Requires title, and ylabel text, and legend text (cell)
    function style_count_plot(title1, xtext, ytext, legtext)
    
    % errorbars can force Y range < 0, so clip them at 0
    set(gca, 'YLim', max(0, get(gca, 'YLim'))); 
    % format the title
    title2 = plot_make_title(t_info);

    % title: prevent special interpretation of _
    title({title2, title1}, 'Interpreter', 'none');
    % axis labels
    xlabel(xtext, 'FontWeight', 'bold')
    ylabel(ytext, 'FontWeight', 'bold')
    % NE corner is good for histograms
    if length(legtext) > 0,
        legend(legtext, 'Interpreter', 'none', 'Location', 'northeast');
    end
    set(gca,'TitleFontSizeMultiplier', 1.1)
    set(gca, 'FontSize', 13);
    grid on;
    end % inner function
    
% Inner function:
%   write the current figure to various files
%   uses the dest_tmpl mfile variable
    function write_plots(dest_name)

    % write the plot out
    if ~isempty(dest_tmpl),
        for ext = ExtList,
            fn_gfx = sprintf(dest_tmpl, dest_name, ext{1});
            fprintf(1, '\tExporting to %s\n', fn_gfx);
            export_fig('-transparent', '-m2', fn_gfx);
        end;
    end
    end

%
% Computation
% 

% histogram domain (in counts)
ct_samp_1 = t_counts{:, 'h_event_count_lo'};
% offsets on bars in the plot -- unused at present
ct_offsets_1 = [0]; % counts units

%% Note: the code below is set up for multiple-line error-bar
%% plots, but we have not used multiple bars.

clf; set(gcf, 'Position', [100 100 850 500]);
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Slews

% just slew count in this plot
names = {'h_event_count_slew'};
names_legend = {'Slews'};
N_plot = length(names);
% indexes to plot
inx = [1:150];

% put the above-selected counts on one plot
for n = 1:N_plot,
    f = names{n};
    f_mean = sprintf('%s_%s', f, 'mean');
    f_std  = sprintf('%s_%s', f, 'std');
    h_eb = plot(ct_samp_1(inx)+ct_offsets_1(n), ...
                       t_counts{inx,f_mean});
    % (error bars just confuse this plot)
    % h_eb{n} = errorbar(ct_samp_1(inx)+ct_offsets_1(n), ...
    %                   t_counts{inx,f_mean}, ...
    %                   t_counts{inx,f_std});
    % style the plot
    set(h_eb, 'LineWidth', 2);
    hold on;
end;

style_count_plot('Mean Event Count: Slews', ...
                 'Number of Events [count]', ...
                 'Frequency [density]', names_legend);
write_plots('event-count-slew');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Detections

clf;

% detections, etc.
names = {'h_event_count_det', 'h_event_count_char', 'h_event_count_detp'};
names_legend = {'Detections', 'Characterizations', 'Detections, Promoted'};
N_plot = length(names);
% indexes to plot (event counts)
inx = [1:500];

% put the above-selected counts on one plot
for n = 1:N_plot,
    f = names{n};
    f_mean = sprintf('%s_%s', f, 'mean');
    f_std  = sprintf('%s_%s', f, 'std');
    % average in bins, standard deviations average in quadrature
    mean_avg = movmean(t_counts{:,f_mean},       3, 'Endpoints', 'fill');
    std_avg = sqrt(movmean(t_counts{:,f_std}.^2, 3, 'Endpoints', 'fill'));
    if false,
        h_eb = errorbar(ct_samp_1(inx)+ct_offsets_1(n), ...
                       mean_avg(inx), ...
                       std_avg(inx));
        set(h_eb, 'LineWidth', 2);
    else,
        plot(ct_samp_1(inx), mean_avg(inx), 'LineWidth', 2);
        % (skip the error bars)
        % hold on;
        % plot(ct_samp_1(inx)+ct_offsets_1(n), ...
        %                mean_avg(inx) + std_avg(inx), 'LineWidth', 1);
    end
    % style the plot
    hold on;
end;

style_count_plot('Mean Event Count: Detections and Characterizations', ...
                 'Number of Events [count]', ...
                 'Frequency [density]', names_legend);
write_plots('event-count-det');

% zoomed
style_count_plot('Mean Event Count: Detections and Characterizations (Zoomed)', ...
                 'Number of Events [count]', ...
                 'Frequency [density]', names_legend);
axis([-Inf 100 0 Inf]);
write_plots('event-count-det-zoom');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Detections - RV precursors

clf;

% detections, etc.
names = {'h_event_count_det_rvplan', 'h_event_count_char_rvplan', 'h_event_count_char'};
names_legend = {'Detections, RV Targets', 'Characterizations, RV Targets', 'Characterizations, Any'};
N_plot = length(names);
% indexes to plot (event counts)
inx = [1:500];

% this guard is needed for old reductions
if ismember('h_event_count_det_rvplan_mean', t_counts.Properties.VariableNames),

    % put the above-selected counts on one plot
    for n = 1:N_plot,
        f = names{n};
        f_mean = sprintf('%s_%s', f, 'mean');
        f_std  = sprintf('%s_%s', f, 'std');
        % average in bins, standard deviations average in quadrature
        mean_avg = movmean(t_counts{:,f_mean},       3, 'Endpoints', 'fill');
        std_avg = sqrt(movmean(t_counts{:,f_std}.^2, 3, 'Endpoints', 'fill'));
        if false,
            h_eb = errorbar(ct_samp_1(inx)+ct_offsets_1(n), ...
                            mean_avg(inx), ...
                            std_avg(inx));
            set(h_eb, 'LineWidth', 2);
        else,
            plot(ct_samp_1(inx), mean_avg(inx), 'LineWidth', 2);
            % (skip the error bars)
            % hold on;
            % plot(ct_samp_1(inx)+ct_offsets_1(n), ...
            %                mean_avg(inx) + std_avg(inx), 'LineWidth', 1);
        end
        % style the plot
        hold on;
    end;

    style_count_plot('Mean Event Count: RV Detections and Characterizations', ...
                     'Number of Events [count]', ...
                     'Frequency [density]', names_legend);
    axis([-Inf 30 0 Inf]);
    write_plots('event-count-det-rv');

end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Earth Characterization Counts

% allow early exit if data not present
if isempty(t_earth_counts),
    return;
end

clf;

% histogram domain (in counts)
% Bins include lower boundary, but exclude upper boundary
% Thus: if "lo" is [0,1,2,...] then first bin is [0,1).
ct_samp_1 = t_earth_counts{:, 'h_earth_char_count_lo'};
% offsets on bars in the plot -- unused at present
ct_offsets_1 = [0 0]; % counts units

% earth char counts appearing in this plot
names = {'h_earth_char_all', 'h_earth_char_strict'};
names_legend = {'Earths Characterized (full/partial, any band)', ...
                'Earths Characterized (full, all bands)'};
N_plot = length(names);

% put the above-selected counts on one plot
for n = 1:N_plot,
    f = names{n};
    f_mean = sprintf('%s_%s', f, 'mean');
    f_std  = sprintf('%s_%s', f, 'std');
    % guard against out-of-date csv files for strict mode
    if ~any(strcmp(f_mean, t_earth_counts.Properties.VariableNames)),
        fprintf('Skipping earth chars (%s): redo make reduce to fix', f_mean);
        continue;
    end
    h_eb = plot(ct_samp_1+ct_offsets_1(n), ...
                       t_earth_counts{:,f_mean});
    % (error bars just confuse this plot)
    % h_eb{n} = errorbar(ct_samp_1+ct_offsets_1(n), ...
    %                   t_earth_counts{:,f_mean}, ...
    %                   t_earth_counts{:,f_std});
    % style the plot
    set(h_eb, 'LineWidth', 2);
    hold on;
end;

% Clip the x-range, ensure 0 gets an explicit bin
axis([-0.5 50 0 Inf]);
style_count_plot('Number of Earths Characterized', ...
                 'Number of Earths [count]', ...
                 'Frequency [density]', names_legend);
write_plots('earth-char-count');

%%% Version 2/a: Showing one mode only, bars

clf;

bar_style = {'FaceColor', [70 130 180]/255}; % SteelBlue

% histogram domain (in counts)
% Bins include lower boundary, but exclude upper boundary
% Thus: if "lo" is [0,1,2,...] then first bin is [0,1).
ct_samp_1 = t_earth_counts{:, 'h_earth_char_count_lo'};
% offsets on bars in the plot -- unused at present
ct_offsets_1 = [0 0]; % counts units

% earth char counts appearing in this plot
names = {'h_earth_char_strict'};
% names_legend = {'Earths Characterized (full, all bands)'};
names_legend = ''; % 1 set of bars -> suppress legend
N_plot = length(names);

% put the above-selected counts on one plot
% *** NB: this bar-plot will fail for N_plot > 1
for n = 1:N_plot,
    f = names{n};
    f_mean = sprintf('%s_%s', f, 'mean');
    f_std  = sprintf('%s_%s', f, 'std');
    % guard against out-of-date csv files for strict mode
    if ~any(strcmp(f_mean, t_earth_counts.Properties.VariableNames)),
        fprintf('Skipping earth chars (%s): redo make reduce to fix', f_mean);
        continue;
    end
    h_eb = bar(ct_samp_1+ct_offsets_1(n), ...
               t_earth_counts{:,f_mean}, bar_style{:});
    % (error bars just confuse this plot)
    % h_eb{n} = errorbar(ct_samp_1+ct_offsets_1(n), ...
    %                   t_earth_counts{:,f_mean}, ...
    %                   t_earth_counts{:,f_std});
    % style the plot
    hold on;
end;

% Clip the x-range
axis([-0.5 50 0 Inf]);
style_count_plot('Number of Earths Characterized', ...
                 'Number of Earths [count]', ...
                 'Frequency [density]', names_legend);
write_plots('earth-char-count-strict');

%%% Version 2/b: Showing one mode only, bars

clf;

% histogram domain (in counts)
% Bins include lower boundary, but exclude upper boundary
% Thus: if "lo" is [0,1,2,...] then first bin is [0,1).
ct_samp_1 = t_earth_counts{:, 'h_earth_char_count_lo'};
% offsets on bars in the plot -- unused at present
ct_offsets_1 = [0 0]; % counts units

% earth char counts appearing in this plot
names = {'h_earth_char_all'};
% names_legend = {'Earths Characterized (full/partial)'};
names_legend = ''; % 1 set of bars -> suppress legend
N_plot = length(names);

% put the above-selected counts on one plot
% *** NB: this bar-plot will fail for N_plot > 1
for n = 1:N_plot,
    f = names{n};
    f_mean = sprintf('%s_%s', f, 'mean');
    f_std  = sprintf('%s_%s', f, 'std');
    % guard against out-of-date csv files for strict mode
    if ~any(strcmp(f_mean, t_earth_counts.Properties.VariableNames)),
        fprintf('Skipping earth chars (%s): redo make reduce to fix', f_mean);
        continue;
    end
    h_eb = bar(ct_samp_1+ct_offsets_1(n), ...
               t_earth_counts{:,f_mean}, bar_style{:});
    % (error bars just confuse this plot)
    % h_eb{n} = errorbar(ct_samp_1+ct_offsets_1(n), ...
    %                   t_earth_counts{:,f_mean}, ...
    %                   t_earth_counts{:,f_std});
    % style the plot
    hold on;
end;

% Clip the x-range
axis([-0.5 50 0 Inf]);
style_count_plot('Number of Earths Characterized', ...
                 'Number of Earths [count]', ...
                 'Frequency [density]', names_legend);
write_plots('earth-char-count-all');


% end of m-file
end
