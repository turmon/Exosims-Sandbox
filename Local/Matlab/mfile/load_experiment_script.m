% load_multi_drms_script -- script for multi-ensemble summaries
%
% * This Matlab script produces some multi-ensemble plots.
%
% Inputs:
%   in_experiment
%   
% Outputs:
%   (to disk)
%
% See also: plot_drms.sh

fprintf('%s: Path info: now running: %s\n', mfilename, mfilename('fullpath'));

% set up output files
if ~exist('dest_tmpl', 'var'),
    dest_tmpl = fullfile(getenv('HOME'), 'expt-det-%s.pdf');
end

% require in_experiment
if ~exist('in_experiment', 'var'),
  error('Require "in_experiment" to be defined');
end;
% check that input exists
if ~exist(in_experiment, 'file'),
    error(sprintf('Could not read input experiment file <%s>', in_experiment));
end

% template for reduced data files: turns a sim-name (and extensions) -> filename
in_tmpl = 'sims/%s/reduce-%s.%s';

expt_info = jsondecode(fileread(in_experiment));

Nscript = length(expt_info.scripts);

% compile experiment into scripts and short tags
%   what we extract: 
script_fns  = cell(Nscript, 1);
script_tags = cell(Nscript, 1);
for s = 1:Nscript,
    script_tags{s} = expt_info.scripts{s}{1};
    script_fns{s} = expt_info.scripts{s}{2};
end;
% extract other info
if isfield(expt_info, 'comment'),
    fprintf('Experiment: %s\n', expt_info.comment);
end;
if isfield(expt_info, 'shortName'),
    expt_tag = sprintf('Experiment: %s', expt_info.shortName);
else
    expt_tag = '';
end;

fprintf('Loading %d different ensembles...\n', Nscript);

% construct by repeated concatenation: inefficient but length is small
%   Ts = struct-array of tables
%   Ts(s).radlum is the radius/luminosity CSV-table for script #s
% NB: if ~exist(...) just allows easy experimentation
if ~exist('Ts', 'var'),
    for s = 1:Nscript,
        t_info   = readtable(sprintf(in_tmpl, script_fns{s}, 'info', 'csv'));
        t_radlum = readtable(sprintf(in_tmpl, script_fns{s}, 'radlum', 'csv'));
        t_earth  = readtable(sprintf(in_tmpl, script_fns{s}, 'earth', 'csv'));
        % holds all the data
        Ts(s) = struct('info', t_info, 'radlum', t_radlum, 'earth', t_earth);
    end;
end;

% re-format Ts into Mradlum, a struct containing matrices with scripts
% across the rows, and planet-bins down the columns.
%  Mradlum.h_RpL_det_main_mean(:,s) is the vectorized histogram for script #s, etc.

% re-format all quantities in the table
% e.g., Qtys = {'h_RpL_det_main_mean', 'h_RpL_det_alt_mean' (etc)};
Qtys = Ts(1).radlum.Properties.VariableNames;

% loop over scripts and quantities
for s = 1:Nscript,
    for qty = Qtys,
        qty1 = qty{1};
        Mradlum.(qty1)(:,s) = Ts(s).radlum{:, qty1};
    end
end

%% Begin graphics

% turn a headless-matlab warning off - a better way might be
% to open a figure window and set "painters" mode in it?
warning('off','MATLAB:prnRenderer:opengl');

figure;

for fig_type = 0:1,

    clf;

    % combo + blue-mode detections...
    %   note: "xdet" includes repeats.  "det" is unique detections only.
    if fig_type == 0,
        tag = 'Unique'; % for text label
        tag_s = 'uniq'; % for filename
        h_RpL_det_mean = Mradlum.h_RpL_det_main_mean + Mradlum.h_RpL_det_alt_mean;
        h_RpL_det_std  = hypot(Mradlum.h_RpL_det_main_std, Mradlum.h_RpL_det_alt_std); % error-bars in quadrature
    else,
        tag = 'Repeat'; % for text label
        tag_s = 'rept'; % for filename
        h_RpL_det_mean = Mradlum.h_RpL_xdet_main_mean + Mradlum.h_RpL_xdet_alt_mean;
        h_RpL_det_std  = hypot(Mradlum.h_RpL_xdet_main_std, Mradlum.h_RpL_xdet_alt_std); % error-bars in quadrature
    end;

    % binnings of planets
    bin_groups = {[1:3], [4:6], [7:9], [10:12], [13:15]};
    bin_names = {'Hot', 'Warm', 'Cold'};
    bin_gnames = {'Rocky', 'Super-Earth', 'Sub-Neptune', 'Sub-Jovian', 'Jovian'};
    bin_gabbr = {'rocky', 'supearth', 'subnep', 'subjup', 'jup'};

    % errorbar plot style
    style_eb = {'LineWidth', 2.5, 'CapSize', 12, 'MarkerSize', 8};

    for b = 1:length(bin_groups),
        clf;
        % the specific exoplanet-bins to plot now
        bin_group1 = bin_groups{b};
        % #planet-bins in the group we're plotting
        Nbin = length(bin_group1); 
        % make the X-pos of each different planet-bin slightly offset for readability
        Xposition = cumsum(ones(Nscript, Nbin)) + ...
            (0.1 * ones(Nscript,1) * [-1:Nbin-2]);
        errorbar(Xposition, ...
                 h_RpL_det_mean(bin_group1, :)', ...
                 h_RpL_det_std(bin_group1, :)', ...
                 'o', style_eb{:});
        legend(bin_names{:});
        % both overall tag and per-plot tag
        title_block = {expt_tag, sprintf('Detections (%s) of %s ExP vs. Ensemble', tag, bin_gnames{b})};
        title(title_block);
        set(gca, 'XTick', [1:Nscript]);
        set(gca, 'XTickLabels', script_tags);
        grid on;
        if length(dest_tmpl) > 0,
            full_tag = sprintf('%s-%s', tag_s, bin_gabbr{b});
            for ext = {'png', 'pdf'},
                fn_gfx = sprintf(dest_tmpl, full_tag, ext{1});
                fprintf('Exporting to "%s"...\n', fn_gfx);
                export_fig('-transparent', '-m2', fn_gfx);
            end;
        end;
    end

end % for fig_type

return

