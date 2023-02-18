% plot_drms_script -- script for ensemble summary plots
%
% * This Matlab script is the top-level Matlab driver for making 
% detection-statistics plots.
% * Inputs are passed by variables in the workspace.
% * The following variables must be defined:
%     in_tmpl:  
%         template-filename for data inputs, containing 2 %s's 
%     dest_tmpl:
%         template-filename for graphical output
%         (the enclosing directory must exist)
% * These variables are optional:
%     mode:
%         a comma-separated list -- currently unused
% 
% * You should probably use the companion shell-script driver (see
% below).  But, if you want to use this manually to debug plots, do:
%   $ matlab -nodesktop
%   >> addpath('Local/Matlab/mfile')
%   >> in_tmpl = 'sims/HabEx_4m_REST_OF_SCRIPT_NAME/reduce-%s.%s';
%   >> t_info = readtable(sprintf(in_tmpl, 'info', 'csv'));
%   >> t_radlum = readtable(sprintf(in_tmpl, 'radlum', 'csv'));
%   >> t_earth = readtable(sprintf(in_tmpl, 'earth', 'csv'));
%   >> mode='';                                                       
%   >> dest_tmpl='';
% and then you can paste lines from this script, for example,
%   >> plot_drm_radlum(dest_tmpl, mode, t_info, t_radlum, t_earth);
%
%
% Inputs:
%   string in_tmpl
%   string dest_tmpl
%   opt string mode = ''
%
% Outputs:
%   (to disk)
%
% See also: plot_drms.sh

fprintf('\n');
fprintf('%s: Path info: now running: %s\n', mfilename, mfilename('fullpath'));
fprintf('%s: Path info: ver is at: %s\n', mfilename, which('ver'));
fprintf('\n');

% require in_tmpl
if ~exist('in_tmpl', 'var'),
  error('Require "in_tmpl" to be defined');
end;
% require dest_tmpl
if ~exist('dest_tmpl', 'var') || length(dest_tmpl) == 0,
  error('Require "dest_tmpl" to be defined and nonempty');
end;
% default mode
if ~exist('mode', 'var'),
  mode = '';
end;

% parse mode
mode_args = strsplit(mode, ';');
mode_param = struct();
mode_param.op = '*'; % default value
for i = 1:length(mode_args),
    if isempty(mode_args{i}),
        continue;
    end;
    mode_arg = strsplit(mode_args{i}, '=');
    if length(mode_arg) ~= 2,
        error('mode argument must be of form FIELD=VALUE separated by ;');
    end;
    if ~isvarname(mode_arg{1}),
        error('illegal fieldname in mode argument');
    end;
    % insert FIELD = VALUE into mode_param
    mode_param.(mode_arg{1}) = mode_arg{2};
end;
    

fprintf('%s: Beginning.\n', mfilename);

% turn a headless-matlab warning off - a better way might be
% to open a figure window and set "painters" mode in it?
warning('off','MATLAB:prnRenderer:opengl');

%%% Useful processing below here

% 1: read the tables - results in t_*

% reduction-info table (basic metadata)
t_info = readtable(sprintf(in_tmpl, 'info', 'csv'));
% time-used table (dets/char/slew/fuel)
t_det_time = readtable(sprintf(in_tmpl, 'times', 'csv'));
% planetary radius/luminosity table
t_radlum = readtable(sprintf(in_tmpl, 'radlum', 'csv'));
% exo-Earth table
t_earth = readtable(sprintf(in_tmpl, 'earth', 'csv'));
% per-star count table
t_star_targ = readtable(sprintf(in_tmpl, 'star-target', 'csv'));
% event duration table
t_events = readtable(sprintf(in_tmpl, 'events', 'csv'));
% event count table
t_event_counts = readtable(sprintf(in_tmpl, 'event-counts', 'csv'));
% earth char list
t_earth_chars = readtable(sprintf(in_tmpl, 'earth-char-list', 'csv'));
% yield-vs-time table
%   yield-vs-time is a superset of old det-vs-time info
t_yield_time = readtable(sprintf(in_tmpl, 'yield-time', 'csv'));
% earth char counts
t_earth_char_count = readtable(sprintf(in_tmpl, 'earth-char-count', 'csv'));
% target promotion table (left if/else because promotion may not always be universal)
if exist(sprintf(in_tmpl, 'promote', 'csv'), 'file'),
    t_promote = readtable(sprintf(in_tmpl, 'promote', 'csv'));
    t_phist   = readtable(sprintf(in_tmpl, 'promote-hist', 'csv'));
else,
    t_promote = [];
    t_phist   = [];
    fprintf('%s: Using null target promotion table: "make reduce" should fix.\n', ...
            mfilename);
end
% visits-vs-time list - can phase out the if/else as processing catches up (2/2023)
if exist(sprintf(in_tmpl, 'visit-time', 'csv'), 'file'),
    t_visit_time = readtable(sprintf(in_tmpl, 'visit-time', 'csv'));
else,
    t_visit_time = [];
    fprintf('%s: Using null visit-vs-time table: "make reduce" should fix.\n', ...
            mfilename);
end

% 2: read the JSON pooled-data file
%    data is available as fields within pool, such as
%        pool.h_xxx_mean, pool.h_xxx_std, ...
%pool = jsondecode(fileread(sprintf(in_tmpl, 'pool', 'json')));

% 3: make the graphics
plot_drm_time_used(      dest_tmpl, mode_param, t_info, t_det_time);
plot_drm_fuel_used(      dest_tmpl, mode_param, t_info, t_det_time);
plot_drm_yield_times(    dest_tmpl, mode_param, t_info, t_yield_time);
plot_drm_visit_times(    dest_tmpl, mode_param, t_info, t_visit_time);
plot_drm_radlum(         dest_tmpl, mode_param, t_info, t_radlum, t_earth);
plot_drm_star_targets(   dest_tmpl, mode_param, t_info, t_star_targ);
plot_drm_events(         dest_tmpl, mode_param, t_info, t_events);
plot_drm_event_counts(   dest_tmpl, mode_param, t_info, t_event_counts, t_earth_char_count);
plot_drm_promote(        dest_tmpl, mode_param, t_info, t_promote, t_phist);
plot_drm_earth_chars(    dest_tmpl, mode_param, t_info, t_earth_chars);
% write the file indicating success
plot_drm_signal_end(dest_tmpl, mode_param, t_info);

%%% Useful processing above here

fprintf('%s: Done.\n', mfilename);

return

