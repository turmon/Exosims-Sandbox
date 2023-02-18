function plot_drm_time_used(dest_tmpl, mode, t_info, t_det_time)
%plot_drm_time_used	plot time-used in a drm-set
% 
% plot_drm_time_used(dest_tmpl, mode, t_info, t_det_time)
% * Time-series plots of time-used by detections, chars, slews, 
% both cumulatively over the mission, and month-by-month.
% * Note: This plot family originally included detection counts (cume,
% uniq, revisit) vs. time - within the t_det_time table.  We later created 
% a newer table (t_yield_time) that has det/char, and allplanets/earth,
% and multiband char info.  That plot family has superseded some of the
% plots that used to be here.
%
% Inputs:
%   string dest_tmpl
%   string mode
%   table t_info
%   table t_det_time
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
   length(strfind(mode.op, 'det_time')) == 0,
    fprintf('Det time plots: skipping, as directed.\n');
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
    function style_det_plot(title1, ytext, legtext)
    
    % errorbars can force Y range < 0, so clip them at 0
    set(gca, 'YLim', max(0, get(gca, 'YLim'))); 
    % format the title
    title2 = plot_make_title(t_info);

    % title: prevent special interpretation of _
    title({title2, title1}, 'Interpreter', 'none');
    % axis labels
    xlabel('Time [day]', 'FontWeight', 'bold')
    ylabel(ytext, 'FontWeight', 'bold')
    legend(legtext, 'Interpreter', 'none', 'Location', 'northwest');
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

% offsets on various error-bars in the plot 
%  (units of days with one complete sample per ~30 days)
t_offsets = [0 10 5];
% sample times
tsamp = t_det_time{:, 'h_det_time_lo'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Slews and Characterization vs. time -- Cumulative

clf; set(gcf, 'Position', [100 100 850 500]);
   
names = {'h_time_det_cume', 'h_time_char_cume', 'h_time_slew_cume'};
names_legend = {'Detection Time', 'Characterization Time', 'Slew Time'};
N_plot = length(names);

% put each of the above detection-times on one plot
h_eb = cell(N_plot, 1);
for n = 1:N_plot,
    f = names{n};
    f_mean = sprintf('%s_%s', f, 'mean');
    f_std  = sprintf('%s_%s', f, 'std');
    h_eb{n} = errorbar(tsamp+t_offsets(n), t_det_time{:,f_mean}, t_det_time{:,f_std});
    % style the plot
    set(h_eb{n}, 'LineWidth', 2);
    hold on;
end;

style_det_plot('Mission Observation Scheduling (Cumulative) vs. Mission Time', ...
               'Cumulative Time Used [day]', names_legend);
write_plots('obstime-cume');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Slews and Characterization vs. time -- Incremental

clf; set(gcf, 'Position', [100 100 850 500]);
   
names = {'h_time_det_incr', 'h_time_char_incr', 'h_time_slew_incr'};
names_legend = {'Detection Time', 'Characterization Time', 'Slew Time'};
N_plot = length(names);

% put each of the above detection-times on one plot
h_eb = cell(N_plot, 1);
for n = 1:N_plot,
    f = names{n};
    f_mean = sprintf('%s_%s', f, 'mean');
    f_std  = sprintf('%s_%s', f, 'std');
    h_eb{n} = errorbar(tsamp+t_offsets(n), t_det_time{:,f_mean}, t_det_time{:,f_std});
    % style the plot
    set(h_eb{n}, 'LineWidth', 1); % Note: thinner
    hold on;
end;

style_det_plot('Mission Observation Scheduling (Incremental) vs. Mission Time', ...
               'Incremental Time Used [days/month]', names_legend);
write_plots('obstime-incr');


return
end % of main function
    