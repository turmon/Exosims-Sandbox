function plot_drm_visit_times(dest_tmpl, mode, t_info, t_visit_time)
%plot_drm_visit_times	plot detection/charaterization visits in a drm-set
% 
% plot_drm_visit_times(dest_tmpl, mode, t_info, t_visit_time)
% * Time-series plots of visits for the purpose of detection or
% characterization - (#visits) plotted vs. mission elapsed time.
% * Note: This plot family is similar to yield, but is *not* yield:
%   t_yield_time -- yield (#planets) versus mission elapsed time. 
%   t_visit_time -- star-visits versus mission elapsed time.
% So, "visits" are counted whether they result in success or not, and 
% visits do not count planets, they count periods-of-integration.
%
% Inputs:
%   string dest_tmpl
%   string mode
%   table t_info
%   table t_visit_time
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
   length(strfind(mode.op, 'visit_time')) == 0,
    fprintf('Visit time plots: skipping, as directed.\n');
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
    function style_visit_plot(title1, ytext, legtext)
    
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
tsamp = t_visit_time{:, 'h_det_time_lo'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot roster

% Needed for each item here: 
%  {fieldnames, plot title, filename}

FigRoster = { ...
    { ...
        {'h_visit_det_visit', 'h_visit_det_uniq', 'h_visit_det_revi'}, ...
        'Detection-Mode Target Visits', 'visit-time-det' ... 
    }, ...
    { ...
        {'h_visit_char_visit', 'h_visit_char_uniq', 'h_visit_char_revi'}, ...
        'Characterization-Mode Target Visits', 'visit-time-char' ... 
    }, ...
             };
    
% iterate over all figures in the roster

for fig_num = 1:length(FigRoster),
    FigInfo = FigRoster{fig_num};
    [names, dtxt, fname] = FigInfo{:};

    % lines in the plot
    N_plot = length(names);

    % skip unless the required fieldnames are present
    % (e.g., no char summaries are in t_visit_time for det-only missions)
    skipping = false;
    for n = 1:N_plot,
        f_mean = sprintf('%s_%s', names{n}, 'mean');
        if ~any(strcmp(fieldnames(t_visit_time), f_mean)),
            skipping = true;
        end
    end;
    if skipping,
        fprintf(sprintf('\tSkipping %s plot (missing field)\n', fname));
        continue;
    end

    %% A: Monthly plot
    clf; set(gcf, 'Position', [100 100 850 500]);
   
    names_legend = {sprintf('All %s', dtxt), ...
                    sprintf('First-time %s', dtxt), ...
                    sprintf('Return %s', dtxt) };

    % put each of the above timeseries on one plot
    for n = 1:N_plot,
        f = names{n};
        f_mean = sprintf('%s_%s', f, 'mean');
        f_std  = sprintf('%s_%s', f, 'std');
        h_eb = errorbar(tsamp+t_offsets(n), t_visit_time{:,f_mean}, t_visit_time{:,f_std});
        % style the plot
        set(h_eb, 'LineWidth', 1); % NB: skinny
        hold on;
    end;

    style_visit_plot(...
        sprintf('Monthly %s vs. Mission Time', dtxt), ...
        sprintf('%s [count/month]', dtxt), ...
        names_legend);
    write_plots([fname '-month']);


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
                        cumsum(t_visit_time{:,f_mean}), sqrt(cumsum(t_visit_time{:,f_std}.^2)));
        % style the plot
        set(h_eb, 'LineWidth', 2);
        hold on;
    end;

    style_visit_plot(...
        sprintf('Cumulative %s vs. Mission Time', dtxt), ...
        sprintf('%s [count]', dtxt), ...
        names_legend);
    write_plots([fname '-cume']);

end % of figure for-loop

return
end % of main function
    