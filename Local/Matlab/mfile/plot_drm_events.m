function plot_drm_events(dest_tmpl, mode, t_info, t_events)
%plot_drm_events	plot event durations in a drm-set
% 
% plot_drm_events(dest_tmpl, mode, t_info, t_events)
% * Plots of event-durations such as characterization integration times,
% slews, etc.
% 
% Inputs:
%   string dest_tmpl
%   string mode
%   table t_info
%   table t_events
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
   length(strfind(mode.op, 'events')) == 0,
    fprintf('Event plots: skipping, as directed.\n');
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
    function style_event_plot(title1, xtext, ytext, legtext)
    
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

% sample times - two resolutions
tsamp_0 = t_events{:, 'h_event_b0_duration_lo'};
tsamp_1 = t_events{:, 'h_event_b1_duration_lo'};
tsamp_2 = t_events{:, 'h_event_b2_duration_lo'};
% offsets on bars in the plot -- unused at present
t_offsets_0 = [0]; % days units
t_offsets_1 = [0]; % days units
t_offsets_2 = t_offsets_1 / 24.0; % days - finer res. (?)

%% Note: the code below is set up for multiple-line error-bar
%% plots, but we have not used multiple bars.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Detection duration
clf; set(gcf, 'Position', [100 100 850 500]);
   
% just detection duration for now
names = {'h_event_det_b1_duration'};
names_legend = {'Detection Time'};
N_plot = length(names);

% put each of the above detection-times on one plot
h_eb = cell(N_plot, 1);
for n = 1:N_plot,
    f = names{n};
    f_mean = sprintf('%s_%s', f, 'mean');
    f_std  = sprintf('%s_%s', f, 'std');
    h_eb{n} = errorbar(tsamp_1+t_offsets_1(n), t_events{:,f_mean}, t_events{:,f_std});
    % style the plot
    set(h_eb{n}, 'LineWidth', 2);
    hold on;
end;

style_event_plot('Mean Integration Time: Detection', ...
                 'Integration Time [day]', ...
                 'Frequency [density]', names_legend);
write_plots('duration-det-b1');

%% same, but hourly binning
clf;
   
% scale factor (days -> hours)
% we're plotting a density, so if we scale X, we have to scale Y 
% to preserve unit area under the density curve
SF = 24; 

names = {'h_event_det_b2_duration'};
names_legend = {'Detection Time'};
N_plot = length(names);

% put each of the above detection-times on one plot
h_eb = cell(N_plot, 1);
for n = 1:N_plot,
    f = names{n};
    f_mean = sprintf('%s_%s', f, 'mean');
    f_std  = sprintf('%s_%s', f, 'std');
    h_eb{n} = errorbar(SF*(tsamp_2+t_offsets_2(n)), t_events{:,f_mean}/SF, t_events{:,f_std}/SF);
    % style the plot
    set(h_eb{n}, 'LineWidth', 2);
    hold on;
end;

style_event_plot('Mean Integration Time: Detection', ...
                 'Integration Time [hour]', ...
                 'Frequency [density]', names_legend);
write_plots('duration-det-b2');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Characterization duration

clf; set(gcf, 'Position', [100 100 850 500]);
   
% char duration
names = {'h_event_char_b1_duration'};
names_legend = {'Characterization Time'};
N_plot = length(names);

% put each of the above time(s) on one plot
h_eb = cell(N_plot, 1);
for n = 1:N_plot,
    f = names{n};
    f_mean = sprintf('%s_%s', f, 'mean');
    f_std  = sprintf('%s_%s', f, 'std');
    h_eb{n} = errorbar(tsamp_1+t_offsets_1(n), t_events{:,f_mean}, t_events{:,f_std});
    % style the plot
    set(h_eb{n}, 'LineWidth', 2);
    hold on;
end;

style_event_plot('Mean Integration Time: Characterization', ...
                 'Integration Time [day]', ...
                 'Frequency [density]', names_legend);
write_plots('duration-char-b1');

%% same, but hourly binning
clf;
   
% scale factor (days -> hours)
SF = 24; 

names = {'h_event_char_b2_duration'};
names_legend = {'Characterization Time'};
N_plot = length(names);

% put each of the above times on one plot
h_eb = cell(N_plot, 1);
for n = 1:N_plot,
    f = names{n};
    f_mean = sprintf('%s_%s', f, 'mean');
    f_std  = sprintf('%s_%s', f, 'std');
    h_eb{n} = errorbar(SF*(tsamp_2+t_offsets_2(n)), t_events{:,f_mean}/SF, t_events{:,f_std}/SF);
    % style the plot
    set(h_eb{n}, 'LineWidth', 2);
    hold on;
end;

style_event_plot('Mean Integration Time: Characterization', ...
                 'Integration Time [hour]', ...
                 'Frequency [density]', names_legend);
write_plots('duration-char-b2');

%% same, but 2day binning
clf;
   
% Scale factor: plot is a density denominated in days,
% but bin-width of x-axis is 2 days.  This makes the
% b1 plot and the b0 plot in the same units (density in days).
SF = 2; 

names = {'h_event_char_b0_duration'};
names_legend = {'Characterization Time'};
N_plot = length(names);

% put each of the above times on one plot
h_eb = cell(N_plot, 1);
for n = 1:N_plot,
    f = names{n};
    f_mean = sprintf('%s_%s', f, 'mean');
    f_std  = sprintf('%s_%s', f, 'std');
    h_eb{n} = errorbar(1*(tsamp_0+t_offsets_0(n)), SF*t_events{:,f_mean}, SF*t_events{:,f_std});
    % style the plot
    set(h_eb{n}, 'LineWidth', 2);
    hold on;
end;

style_event_plot('Mean Integration Time: Characterization', ...
                 'Integration Time [day]', ...
                 'Frequency [density]', names_legend);
write_plots('duration-char-b0');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Slew duration

clf; set(gcf, 'Position', [100 100 850 500]);
   
% slew duration
names = {'h_event_slew_b1_duration'};
names_legend = {'Slew Time'};
N_plot = length(names);

% put each of the above time(s) on one plot
h_eb = cell(N_plot, 1);
for n = 1:N_plot,
    f = names{n};
    f_mean = sprintf('%s_%s', f, 'mean');
    f_std  = sprintf('%s_%s', f, 'std');
    h_eb{n} = errorbar(tsamp_1+t_offsets_1(n), t_events{:,f_mean}, t_events{:,f_std});
    % style the plot
    set(h_eb{n}, 'LineWidth', 2);
    hold on;
end;

style_event_plot('Mean Times: Slew', ...
                 'Slew Time [day]', ...
                 'Frequency [density]', names_legend);
write_plots('duration-slew-b1');

%% same, but two-daily binning
clf;
   
% slew duration
names = {'h_event_slew_b0_duration'};
names_legend = {'Slew Time'};
N_plot = length(names);

% Scale factor: plot is a density denominated in days,
% but bin-width of x-axis is 2 days.  This makes the
% b1 plot and the b0 plot in the same units (density in days).
SF = 2;

% put each of the above time(s) on one plot
h_eb = cell(N_plot, 1);
for n = 1:N_plot,
    f = names{n};
    f_mean = sprintf('%s_%s', f, 'mean');
    f_std  = sprintf('%s_%s', f, 'std');
    h_eb{n} = errorbar(tsamp_0+t_offsets_0(n), SF*t_events{:,f_mean}, SF*t_events{:,f_std});
    % style the plot
    set(h_eb{n}, 'LineWidth', 2);
    hold on;
end;

style_event_plot('Mean Times: Slew', ...
                 'Slew Time [day]', ...
                 'Frequency [density]', names_legend);
write_plots('duration-slew-b0');

return

% end of m-file
end
