function plot_drm_fuel_use(dest_tmpl, mode, t_info, t_fuel)
%plot_drm_fuel_use	plot fuel use versus time in a drm-set
% 
% plot_drm_fuel_use(dest_tmpl, mode, t_info, t_fuel)
% * Time-series plots of fuel use.  The "t_fuel" input is the same as
% "t_det_time" elsewhere -- it contains more than just fuel.
% 
% Inputs:
%   string dest_tmpl
%   string mode
%   table t_info
%   table t_fuel
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

% file extensions to write
%ExtList = {'png', 'pdf'};
ExtList = {'png'};

%
% Computation
% 

if ~any(strcmp('h_time_fuel_slew_mean', fieldnames(t_fuel))),
    fprintf('No fuel use in the given table, skipping.\n');
    return;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Fuel (in kg) plot
%%

% best to set up the size (aspect ratio is last two numbers)
clf;
set(gcf, 'Position', [100 100 850 500]);
   
% what to plot
names = {'h_time_fuel_slew', 'h_time_fuel_keep', 'h_time_fuel_all'};
names_legend = {'Slew', 'Stationkeeping', 'Total'};
t_offsets = [0 10 5];
N_plot = length(names);

% sample times
tsamp = t_fuel{:, 'h_det_time_lo'};

% put each of the above on one plot
h_eb = cell(N_plot, 1);
top = 0;
for n = 1:N_plot,
    f = names{n};
    f_mean = sprintf('%s_%s', f, 'mean');
    f_std  = sprintf('%s_%s', f, 'std');
    h_eb{n} = errorbar(tsamp + t_offsets(n), t_fuel{:,f_mean}, t_fuel{:,f_std});
    % style the plot
    set(h_eb{n}, 'LineWidth', 2);
    hold on;
    % get some control over axis range
    top = max(top, max(t_fuel{:,f_mean}));
end;

%% Plot/axis styles

title1 = 'Cumulative Fuel Use vs. Time';
title2 = plot_make_title(t_info);

% control axis y-range:
%   negative not allowed
%   top of range can't be too big
%   (if top guards against zero fuel use for empty DRMs)
if top > 0,
    a_range = axis;
    a_range(3) = max(0, a_range(3)); % fuel is >= 0
    a_range(4) = 1.2*top; % some factor above the largest mean
    axis(a_range);
end;

% title: prevent special interpretation of _
title({title2, title1}, 'Interpreter', 'none');
xlabel('Time [days]', 'FontWeight', 'bold')
ylabel('Fuel Used [kg]', 'FontWeight', 'bold')
legend(names_legend, 'Location', 'northwest', 'Interpreter', 'none')
set(gca, 'TitleFontSizeMultiplier', 1.1)
set(gca, 'FontSize', 13);
grid on;


%% write the plot out
if ~isempty(dest_tmpl),
    for ext = ExtList,
        fn_gfx = sprintf(dest_tmpl, 'fuel', ext{1});
        fprintf(1, '\tExporting to %s\n', fn_gfx);
        export_fig('-transparent', '-m2', fn_gfx);
    end;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Delta-v plot
%%

if ~any(strcmp('h_time_delta_v_slew_cume_mean', fieldnames(t_fuel))),
    fprintf('No delta-v in the given table, skipping.  Redo reduce to fix.\n');
    return;
end;

% best to set up the size (aspect ratio is last two numbers)
clf;
set(gcf, 'Position', [100 100 850 500]);
   
% what to plot
% (old)
%names = {'h_time_delta_v_slew_cume', 'h_time_delta_v_obs_cume', 'h_time_delta_v_all_cume'};
%names_legend = {'Slew', 'Observing', 'Total'};
% "total" makes no sense for starshades (two separate objects moving),
% so do not display the sum (.._all_...) as a separate line.
names = {'h_time_delta_v_slew_cume', 'h_time_delta_v_obs_cume'};
names_legend = {'Slew', 'Observing'};
t_offsets = [0 10 5];
N_plot = length(names);

% sample times
tsamp = t_fuel{:, 'h_det_time_lo'};

% put each of the above on one plot
h_eb = cell(N_plot, 1);
top = 0;
for n = 1:N_plot,
    f = names{n};
    f_mean = sprintf('%s_%s', f, 'mean');
    f_std  = sprintf('%s_%s', f, 'std');
    h_eb{n} = errorbar(tsamp + t_offsets(n), t_fuel{:,f_mean}, t_fuel{:,f_std});
    % style the plot
    set(h_eb{n}, 'LineWidth', 2);
    hold on;
    % get some control over axis range
    top = max(top, max(t_fuel{:,f_mean}));
end;

%% Plot/axis styles

title1 = 'Cumulative Delta-V vs. Time';
title2 = plot_make_title(t_info);

% control axis y-range:
%   negative not allowed
%   top of range can't be too big
%   (if top guards against zero fuel use for empty DRMs)
if top > 0,
    a_range = axis;
    a_range(3) = max(0, a_range(3)); % delta-v is >= 0
    a_range(4) = 1.2*top; % some factor above the largest mean
    axis(a_range);
end;

% title: prevent special interpretation of _
title({title2, title1}, 'Interpreter', 'none');
xlabel('Time [days]', 'FontWeight', 'bold')
ylabel('Delta V [m/s]', 'FontWeight', 'bold')
legend(names_legend, 'Location', 'northwest', 'Interpreter', 'none')
set(gca, 'TitleFontSizeMultiplier', 1.1)
set(gca, 'FontSize', 13);
grid on;


%% write the plot out
if ~isempty(dest_tmpl),
    for ext = ExtList,
        fn_gfx = sprintf(dest_tmpl, 'delta-v', ext{1});
        fprintf(1, '\tExporting to %s\n', fn_gfx);
        export_fig('-transparent', '-m2', fn_gfx);
    end;
end

return
