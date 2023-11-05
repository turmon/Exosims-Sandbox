function plot_drm_yield_times(dest_tmpl, mode, t_info, t_yield_time)
%plot_drm_yield_times	plot detection/charaterization times in a drm-set
% 
% plot_drm_yield_times(dest_tmpl, mode, t_info, t_yield_time)
% * Time-series plots of detections, characterizations - plotted vs. 
% mission elapsed time.
% * Note: This plot family originally included only detections (cume,
% uniq, revisit) vs. time - within the t_det_time table.  We later created 
% a newer table (t_yield_time) that has det/char, and allplanets/earth,
% and multiband char info.  This latter table is the one we're using.
%
% Inputs:
%   string dest_tmpl
%   string mode
%   table t_info
%   table t_yield_time
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
   length(strfind(mode.op, 'yield_time')) == 0,
    fprintf('Yield time plots: skipping, as directed.\n');
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
    function style_yield_plot(title1, ytext, legtext)
    
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
tsamp = t_yield_time{:, 'h_det_time_lo'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Detections vs. time

% Needed for each item here: 
%  {fieldnames, plot title, filename}

% NOTE: 
%  * The third plot entry is *only full chars*.  
%  * The last two plot entries *include full and partial chars*.  
% I.e., both full and partial chars are included in the latter two.  
% This extra detail in the plot title was removed to
% make the plots in the SDET report less overwhelming.
FigRoster = { ...
    { ...
        {'h_time_det_allplan_cume', 'h_time_det_allplan_uniq', 'h_time_det_allplan_revi'}, ...
        'Detections [All]', 'time-det-allplan' ... 
    }, ...
    { ...
        {'h_time_det_earth_cume', 'h_time_det_earth_uniq', 'h_time_det_earth_revi'}, ...
        'Detections [Earths]', 'time-det-earth' ... 
    }, ...
    { ...
        {'h_time_char_full_allplan_cume_union', ...
         'h_time_char_full_allplan_uniq_union', 'h_time_char_full_allplan_revi_union'}, ...
        'Characterizations [All Planets]', 'time-char-allplan-full' ...
    }, ...
    { ...
        {'h_time_char_part_allplan_cume_union', ...
         'h_time_char_part_allplan_uniq_union', 'h_time_char_part_allplan_revi_union'}, ...
        'Characterizations [All Planets]', 'time-char-allplan-part' ...
    }, ...
    { ...
        {'h_time_char_part_earth_cume_union', ...
         'h_time_char_part_earth_uniq_union', 'h_time_char_part_earth_revi_union'}, ...
        'Characterizations [Earths]', 'time-char-earth-part' ...
    } ...
             };
    
% 2023-11: skip the montly a/k/a incremental plots
DO_INCREMENTAL_PLOT = false;

% iterate over all figures in the roster

for fig_num = 1:length(FigRoster),
    FigInfo = FigRoster{fig_num};
    [names, dtxt, fname] = FigInfo{:};

    % lines in the plot
    N_plot = length(names);

    % skip unless the required fieldnames are present
    % (e.g., no char summaries are in t_yield_time for det-only missions)
    skipping = false;
    for n = 1:N_plot,
        f_mean = sprintf('%s_%s', names{n}, 'mean');
        if ~any(strcmp(fieldnames(t_yield_time), f_mean)),
            skipping = true;
        end
    end;
    if skipping,
        fprintf(sprintf('\tSkipping %s plot (missing field)\n', fname));
        continue;
    end

    % names = {'h_time_det_allplan_cume', 'h_time_det_allplan_uniq', 'h_time_det_allplan_revi'};
    names_legend = {sprintf('All %s', dtxt), ...
                    sprintf('Unique %s', dtxt), ...
                    sprintf('Revisit %s', dtxt) };

    %% A: Monthly plot
    if DO_INCREMENTAL_PLOT,
        clf; set(gcf, 'Position', [100 100 850 500]);
        
        % put each of the above timeseries on one plot
        for n = 1:N_plot,
            f = names{n};
            f_mean = sprintf('%s_%s', f, 'mean');
            f_std  = sprintf('%s_%s', f, 'std');
            h_eb = errorbar(tsamp+t_offsets(n), t_yield_time{:,f_mean}, t_yield_time{:,f_std});
            % style the plot
            set(h_eb, 'LineWidth', 1); % NB: skinny
            hold on;
        end;

        style_yield_plot(...
            sprintf('Monthly %s vs. Mission Time', dtxt), ...
            sprintf('%s [count/month]', dtxt), ...
            names_legend);
        write_plots([fname '-month']);
    end;


    %% B: Cumulative plot
    % Note: we derive the error bars here

    clf; set(gcf, 'Position', [100 100 850 500]);

    % put each of the above timeseries on one plot
    for n = 1:N_plot,
        f = names{n};
        f_mean = sprintf('%s_%s', f, 'mean');
        f_std  = sprintf('%s_%s', f, 'std');
        % cumulative mean, plus cumulative std error (adding in quadrature,
        % i.e., sqrt-sum-of-squares)
        h_eb = errorbar(tsamp+t_offsets(n), ...
                        cumsum(t_yield_time{:,f_mean}), sqrt(cumsum(t_yield_time{:,f_std}.^2)));
        % style the plot
        set(h_eb, 'LineWidth', 2);
        hold on;
    end;

    style_yield_plot(...
        sprintf('Cumulative %s vs. Mission Time', dtxt), ...
        sprintf('%s [count]', dtxt), ...
        names_legend);
    write_plots([fname '-cume']);

end % of figure for-loop

return
end % of main function
    
