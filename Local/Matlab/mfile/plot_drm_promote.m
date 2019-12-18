function plot_drm_promote(dest_tmpl, mode, t_info, t_promote, t_phist)
%plot_drm_promote	plot target promotion counts in a drm-set
% 
% plot_drm_promote(dest_tmpl, mode, t_info, t_promote)
% * Time-series plots of target promotion criteria
% * Histograms of planet-counts meeting criteria
% 
% Inputs:
%   string dest_tmpl
%   string mode
%   table t_info
%   table t_promote
%   table t_phist
% 
% Outputs:
%   (to disk)
% 
% See Also:  

% Skip plots if the input was []
if ismatrix(t_promote) && isempty(t_promote),
   fprintf('Promotion plots: No data.  Skipping.\n');
   return;
end

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
   length(strfind(mode.op, 'promote')) == 0,
    fprintf('Promotion plots: skipping, as directed.\n');
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
    function style_promote_plot(title1, ytext, legtext)
    
    % errorbars can force Y range < 0, so clip them at 0
    set(gca, 'YLim', max(0, get(gca, 'YLim'))); 
    % format the title
    title2 = plot_make_title(t_info);

    % title: prevent special interpretation of _
    title({title2, title1}, 'Interpreter', 'none');
    % axis labels
    xlabel('Detection Integration Time [day]', 'FontWeight', 'bold')
    ylabel(ytext, 'FontWeight', 'bold')
    legend(legtext, 'Interpreter', 'none', 'Location', 'northwest');
    set(gca,'TitleFontSizeMultiplier', 1.1)
    set(gca, 'FontSize', 13);
    grid on;
    end % inner function

% Inner function:
%   Set up plot/axis styles, title, axis labels.
%   Uses global mfile variables, such as t_info
%   Requires title and legend text (cell)
    function style_phist_plot(title1, legtext)
    
    % errorbars can force Y range < 0, so clip them at 0
    set(gca, 'YLim', max(0, get(gca, 'YLim'))); 
    % format the title
    title2 = plot_make_title(t_info);

    % title: prevent special interpretation of _
    title({title2, title1}, 'Interpreter', 'none');
    % axis labels
    xlabel('Planet Count [count]', 'FontWeight', 'bold')
    ylabel('Density', 'FontWeight', 'bold')
    legend(legtext, 'Interpreter', 'none', 'Location', 'northwest');
    set(gca,'TitleFontSizeMultiplier', 1.1)
    set(gca, 'FontSize', 13);
    set(get(gca, 'Children'), 'MarkerSize', 8, 'LineWidth', 2);
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
%  (fourth offset is rarely used)
t_offsets = [0 10 5 -5];
% sample times
tsamp = t_promote{:, 'h_promo_time_lo'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Promotions vs. time

clf; set(gcf, 'Position', [100 100 850 500]);

%% A: Cumulative AllPlan plot
clf;
names = {'h_promo_count_allplan_cume', 'h_promo_span_allplan_cume', 'h_promo_promo_allplan_cume'};
names_legend = {'Obs. Count', 'Obs. Span', 'Promotion'};
N_plot = length(names);

% put each of the above detection-times on one plot
for n = 1:N_plot,
    f = names{n};
    f_mean = sprintf('%s_%s', f, 'q50');
    f_bar1 = sprintf('%s_%s', f, 'q25');
    f_bar2 = sprintf('%s_%s', f, 'q75');
    h_eb = errorbar(tsamp+t_offsets(n), ...
                    t_promote{:,f_mean}, ...
                    t_promote{:,f_mean} - t_promote{:,f_bar1}, ...
                    t_promote{:,f_bar2} - t_promote{:,f_mean});
    % style the plot
    set(h_eb, 'LineWidth', 1); % NB: skinny
    hold on;
end;

style_promote_plot('All Planets: Cumulative Promotions vs. Detector Time', ...
               'Targets Passing Criterion [count]', names_legend);
write_plots('promote-allplan-cume');

%% B: Cumulative Hab Zone plot
clf;
names = {'h_promo_count_hzone_cume', 'h_promo_span_hzone_cume', 'h_promo_promo_hzone_cume'};
names_legend = {'Obs. Count', 'Obs. Span', 'Promotion'};
N_plot = length(names);

% put each of the above detection-times on one plot
for n = 1:N_plot,
    f = names{n};
    f_mean = sprintf('%s_%s', f, 'q50');
    f_bar1 = sprintf('%s_%s', f, 'q25');
    f_bar2 = sprintf('%s_%s', f, 'q75');
    h_eb = errorbar(tsamp+t_offsets(n), ...
                    t_promote{:,f_mean}, ...
                    t_promote{:,f_mean} - t_promote{:,f_bar1}, ...
                    t_promote{:,f_bar2} - t_promote{:,f_mean});
    % style the plot
    set(h_eb, 'LineWidth', 1); % NB: skinny
    hold on;
end;

style_promote_plot('Habitable Zone: Cumulative Promotions vs. Detector Time', ...
               'Targets Passing Criterion [count]', names_legend);
write_plots('promote-hzone-cume');

%% C: Cumulative Earth plot
clf;
names = {'h_promo_count_earth_cume', 'h_promo_span_earth_cume', 'h_promo_promo_earth_cume'};
names_legend = {'Obs. Count', 'Obs. Span', 'Promotion'};
N_plot = length(names);

% put each of the above detection-times on one plot
for n = 1:N_plot,
    f = names{n};
    f_mean = sprintf('%s_%s', f, 'q50');
    f_bar1 = sprintf('%s_%s', f, 'q25');
    f_bar2 = sprintf('%s_%s', f, 'q75');
    h_eb = errorbar(tsamp+t_offsets(n), ...
                    t_promote{:,f_mean}, ...
                    t_promote{:,f_mean} - t_promote{:,f_bar1}, ...
                    t_promote{:,f_bar2} - t_promote{:,f_mean});
    % style the plot
    set(h_eb, 'LineWidth', 1); % NB: skinny
    hold on;
end;

style_promote_plot('Earthlike Planets: Cumulative Promotions vs. Detector Time', ...
               'Targets Passing Criterion [count]', names_legend);
write_plots('promote-earth-cume');


%% D: Cumulative Star plot
clf;
names = {'h_promo_count_star_cume', 'h_promo_span_star_cume', 'h_promo_promo_star_cume'};
names_legend = {'Obs. Count', 'Obs. Span', 'Promotion'};
N_plot = length(names);

% put each of the above detection-times on one plot
for n = 1:N_plot,
    f = names{n};
    f_mean = sprintf('%s_%s', f, 'q50');
    f_bar1 = sprintf('%s_%s', f, 'q25');
    f_bar2 = sprintf('%s_%s', f, 'q75');
    h_eb = errorbar(tsamp+t_offsets(n), ...
                    t_promote{:,f_mean}, ...
                    t_promote{:,f_mean} - t_promote{:,f_bar1}, ...
                    t_promote{:,f_bar2} - t_promote{:,f_mean});
    % style the plot
    set(h_eb, 'LineWidth', 1); % NB: skinny
    hold on;
end;

style_promote_plot('Counted By Star: Earthlike Planets: Cumulative Promotions vs. Detector Time', ...
               'Targets Passing Criterion [count]', names_legend);
write_plots('promote-star-cume');

%% E: Cumulative Span options plot
clf;
names = {'h_promo_spanPlan_star_cume', 'h_promo_spanHZ0_star_cume', 'h_promo_spanHZ1_star_cume', 'h_promo_spanEarth_star_cume'};
names_legend = {'Span: Any Planet Period', 'Span: Inner HZ Period', 'Span: Outer HZ Period', 'Span: Earthlike Period'};
N_plot = length(names);

% put each of the above detection-times on one plot
for n = 1:N_plot,
    f = names{n};
    f_mean = sprintf('%s_%s', f, 'q50');
    f_bar1 = sprintf('%s_%s', f, 'q25');
    f_bar2 = sprintf('%s_%s', f, 'q75');
    h_eb = errorbar(tsamp+t_offsets(n), ...
                    t_promote{:,f_mean}, ...
                    t_promote{:,f_mean} - t_promote{:,f_bar1}, ...
                    t_promote{:,f_bar2} - t_promote{:,f_mean});
    % style the plot
    set(h_eb, 'LineWidth', 1); % NB: skinny
    hold on;
end;

style_promote_plot('Counted By Star: Cumulative Observation Span vs. Detector Time', ...
               'Targets Passing Span Criterion [count]', names_legend);
write_plots('promote-star-span-cume');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Promotion-count histograms (function of planet count)

x_values = t_phist{:,'h_phist_count_lo'};
names_legend = {'Obs. Count', 'Obs. Span', 'Promotion'};
% line-styles to distinguish often-overlapping lines
LSct = 'o-'; % ct = count
LSsp = 'o-'; % sp = span
LSpr = 'x-'; % pr = promo
LSx  = '+-'; % extra
% small offset to one line
del = 0.00;

%% A: HZ Plot

clf;

title_tmpl = '%s Planets: Promotion Tiers [%s Mission Time]';

plot(x_values, t_phist{:,'h_phist_t1_count_hzone_mean'}+del, LSct, ...
     x_values, t_phist{:,'h_phist_t1_span_hzone_mean'}, LSsp, ...
     x_values, t_phist{:,'h_phist_t1_promo_hzone_mean'}, LSpr);
style_phist_plot(sprintf(title_tmpl, 'Habitable Zone', '2 Years'), names_legend);
write_plots('phist-hzone-2year');

clf;

plot(x_values, t_phist{:,'h_phist_t2_count_hzone_mean'}+del, LSct, ...
     x_values, t_phist{:,'h_phist_t2_span_hzone_mean'}, LSsp, ...
     x_values, t_phist{:,'h_phist_t2_promo_hzone_mean'}, LSpr);
style_phist_plot(sprintf(title_tmpl, 'Habitable Zone', '3 Years'), names_legend);
write_plots('phist-hzone-3year');


%% B: Earth plot

clf;

plot(x_values, t_phist{:,'h_phist_t1_count_earth_mean'}+del, LSct, ...
     x_values, t_phist{:,'h_phist_t1_span_earth_mean'}, LSsp, ...
     x_values, t_phist{:,'h_phist_t1_promo_earth_mean'}, LSpr);
style_phist_plot(sprintf(title_tmpl, 'Earthlike', '2 Years'), names_legend);
write_plots('phist-earth-2year');

clf;

plot(x_values, t_phist{:,'h_phist_t2_count_earth_mean'}+del, LSct, ...
     x_values, t_phist{:,'h_phist_t2_span_earth_mean'}, LSsp, ...
     x_values, t_phist{:,'h_phist_t2_promo_earth_mean'}, LSpr);
style_phist_plot(sprintf(title_tmpl, 'Earthlike', '3 Years'), names_legend);
write_plots('phist-earth-3year');

%% C: Counting-by-Star plot

clf;

plot(x_values, t_phist{:,'h_phist_t1_count_star_mean'}+del, LSct, ...
     x_values, t_phist{:,'h_phist_t1_span_star_mean'}, LSsp, ...
     x_values, t_phist{:,'h_phist_t1_promo_star_mean'}, LSpr);
style_phist_plot(sprintf(title_tmpl, 'By Star, with Earthlike', '2 Years'), names_legend);
write_plots('phist-star-2year');

clf;

plot(x_values, t_phist{:,'h_phist_t2_count_star_mean'}+del, LSct, ...
     x_values, t_phist{:,'h_phist_t2_span_star_mean'}, LSsp, ...
     x_values, t_phist{:,'h_phist_t2_promo_star_mean'}, LSpr);
style_phist_plot(sprintf(title_tmpl, 'By Star, with Earthlike', '3 Years'), names_legend);
write_plots('phist-star-3year');

%% D: Span plot

clf;

title_tmpl_span = '%s Planets: Various Span Critera Only [%s Mission Time]';
names_legend_span = {'Obs. Span > T/2 (Any Actual Planet)', 'Obs. Span > T/2 (Inner HZ, Hypothet.)', ...
                    'Obs. Span > T/2 (Outer HZ, Hypothet.)', 'Obs. Span > T/2 (Actual Earthlike)'};

plot(x_values, t_phist{:,'h_phist_t1_spanPlan_star_mean'}+del, LSct, ...
     x_values, t_phist{:,'h_phist_t1_spanHZ0_star_mean'}, LSsp, ...
     x_values, t_phist{:,'h_phist_t1_spanHZ1_star_mean'}, LSpr, ...
     x_values, t_phist{:,'h_phist_t1_spanEarth_star_mean'}, LSx);
style_phist_plot(sprintf(title_tmpl_span, 'By Star, with Earthlike', '2 Years'), names_legend_span);
write_plots('phist-star-span-2year');

clf;

plot(x_values, t_phist{:,'h_phist_t2_spanPlan_star_mean'}+del, LSct, ...
     x_values, t_phist{:,'h_phist_t2_spanHZ0_star_mean'}, LSsp, ...
     x_values, t_phist{:,'h_phist_t2_spanHZ1_star_mean'}, LSpr, ...
     x_values, t_phist{:,'h_phist_t2_spanEarth_star_mean'}, LSx);
style_phist_plot(sprintf(title_tmpl_span, 'By Star, with Earthlike', '3 Years'), names_legend_span);
write_plots('phist-star-span-3year');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return
end % of main function
    
