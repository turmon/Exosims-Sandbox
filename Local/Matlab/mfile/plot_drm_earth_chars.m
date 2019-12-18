function plot_drm_earth_chars(dest_tmpl, mode, t_info, t_chars)
%plot_drm_earth_chars	plot attempted earth characterizations in a drm-set
% 
% plot_drm_earth_chars(dest_tmpl, mode, t_info, t_chars)
% * Plots of attempted characterizations, e.g., fails/successes in
% WA/dMag coordinates.
% 
% Inputs:
%   string dest_tmpl
%   string mode
%   table t_info
%   table t_chars
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
   length(strfind(mode.op, 'earth-chars')) == 0,
    fprintf('Earth char attempt plots: skipping, as directed.\n');
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
    function style_wa_dmag_plot(title1, bartext)
    
    % errorbars can force Y range < 0, so clip them at 0
    set(gca, 'YLim', max(0, get(gca, 'YLim'))); 
    % format the title
    title2 = plot_make_title(t_info);

    % title: prevent special interpretation of _
    title({title2, title1}, 'Interpreter', 'none');
    % axis labels
    xlabel('Working Angle (mas)', 'FontWeight', 'bold')
    ylabel('Delta Magnitude', 'FontWeight', 'bold')
    % SE corner is good for wa/dmag
    if length('legtext') < 0,
        legend('legtext', 'Interpreter', 'none', 'Location', 'southeast');
    end
    set(gca,'TitleFontSizeMultiplier', 1.1)
    set(gca, 'FontSize', 13);
    hc = colorbar;
    if length(bartext) > 0,
        ylabel(hc, bartext, 'FontWeight', 'bold');
    end;
    grid on;
    % "easy" rectangle: 
    %   upper-left corner at: (x=70 mas, y=26 mag)
    %   lower-left corner at  (x=70 mas, y=ax(3)) 
    ax = axis;
    pR = rectangle('Position', [70 ax(3) ax(2)-70 26-ax(3)]);
    set(pR,'EdgeColor',[1 0 0],'LineStyle','--','LineWidth',1);
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
            % -m is magnification factor, increased above 2 for these dense plots
            % (note: -m3 gives a warning message about a very large
            % intermediate image)
            export_fig('-transparent', '-m2.5', fn_gfx);
        end;
    end
    end

%
% Computation
% 

% for missions without chars, the table will be basically empty.
% skip such tables, there is nothing to do.
if ~any(strcmp(fieldnames(t_chars), 'is_success')),
    fprintf('\tNo char info, skipping earth chars plots\n');
    return;
end;

% booleans
ok    =  t_chars{:, 'is_success'} >  0;
ok2   = (t_chars{:, 'is_success'} == 0) & (t_chars{:, 'n_success'} > 0);
fail  = ~(ok | ok2);
deep  = t_chars{:, 'is_deep'}  > 0;
promo = (t_chars{:, 'is_promo'} > 0) & (~deep);
% real numbers
ec_wa    = t_chars{:, 'WA'};
ec_dmag  = t_chars{:, 'dMag'};
ec_snr   = t_chars{:, 'char_SNR'};
ec_phi   = t_chars{:, 'phi'};
% not all generated CSVs have this field (3/2019)
if any(strcmp('MV', fieldnames(t_chars))),
    ec_vmag  = t_chars{:, 'MV'};
else,
    ec_vmag  = zeros(size(ec_wa));
end
    
% derived
easy = (ec_dmag < 26) & (ec_wa > 70);

%% Table of wa/dmag plots to make.
% Format: textual full name, abbreviated name, index variable
dprops = {'SizeData', 300};
pprops = {'SizeData', 100, 'MarkerEdgeAlpha', 0.5};

PlotRoster = { ...
    {'Deep Dive', 'deep',  deep,  dprops}, ...
    {'Promoted',  'promo', promo, pprops} ...
             };

% create new figures for each plot?
new_figs = true;

for pnum = 1:length(PlotRoster),
    full_name = PlotRoster{pnum}{1};
    abbr_name = PlotRoster{pnum}{2};
    Cinx      = PlotRoster{pnum}{3};
    gprops    = PlotRoster{pnum}{4};
    
    % WA/dMag, shading = char SNR
    if new_figs; figure(); end; clf;
    set(gcf, 'Position', [100 100 850 500]);
    scatter(ec_wa(ok), ec_dmag(ok), 15, [0 0 0]+0.5, '+');
    hold on
    % scatter(ec_wa(~ok & promo), ec_dmag(~ok & promo), 20, ec_snr(~ok&promo),'s');
    pS = scatter(ec_wa(fail & Cinx), ec_dmag(fail & Cinx), 400, ...
                 ec_snr(fail & Cinx), '.', gprops{:});
    pS2= scatter(ec_wa(ok2 & Cinx), ec_dmag(ok2 & Cinx), 15, ...
                 ec_snr(ok2 & Cinx), 'x');
    ax = axis;
    txt_x = ax(1:2) * [0.5; 0.5]; % x midline for text annotations
    txt_y = ax(3:4) * [0.7; 0.3]; % y midline for text annotations

    style_wa_dmag_plot( ...
        sprintf('%s: Earth Characterizations vs. WA and dMag, shaded by Char. SNR', full_name), ...
        'Char SNR');
    legend({sprintf('Successful Chars (%d)', sum(ok)), ...
            sprintf('Failed %s Chars (%d)', full_name, sum(fail & Cinx)), ...
            sprintf('Split Multi-Earth %s Chars (%d)', full_name, sum(ok2  & Cinx))}, ...
           'Location', 'northeast');
    % easy chars sub-legend
    text(txt_x, txt_y, { ...
        '"Easy" Characterizations within red box:', ...
        sprintf('  %d Successful', sum(easy & ok)), ...
        sprintf('  %d Failed (counting %s only)', sum(easy & Cinx & fail), full_name)}, ...
              'FontSize', 13);
    write_plots(sprintf('earth-char-wa-dmag-%s-snr', abbr_name));

    % WA/dMag, shading = log10(Phi)
    if new_figs; figure(); end; clf;
    set(gcf, 'Position', [100 100 850 500]);
    scatter(ec_wa(ok), ec_dmag(ok), 15, [0 0 0]+0.5, '+');
    hold on
    % scatter(ec_wa(~ok & promo), ec_dmag(~ok & promo), 20, log10(ec_phi(~ok&promo)),'s');
    pS = scatter(ec_wa(~ok & Cinx), ec_dmag(~ok & Cinx), 400, ...
            log10(ec_phi(~ok & Cinx)), '.', gprops{:});
    ax = axis;

    style_wa_dmag_plot(...
        sprintf('%s: Earth Characterizations vs. WA and dMag, shaded by Log Phi', full_name),...
        'log_{10}(\Phi) : truncated at -2');
    legend({sprintf('Successful Chars (%d)', sum(ok)), ...
            sprintf('Failed %s Chars (%d)', full_name, sum(~ok & Cinx))}, ...
           'Location', 'northeast');
    % easy chars sub-legend
    text(txt_x, txt_y, { ...
        '"Easy" Characterizations within red box:', ...
        sprintf('  %d Successful', sum(easy & ok)), ...
        sprintf('  %d Failed (counting %s only)', sum(easy & Cinx & ~ok), full_name)}, ...
              'FontSize', 13);
    % -2 on log10 scale = 0.01 = 5 magnitudes
    caxis([-2 0]);
    write_plots(sprintf('earth-char-wa-dmag-%s-phi', abbr_name));

    % WA/dMag, shading = Vmag
    if new_figs; figure(); end; clf;
    set(gcf, 'Position', [100 100 850 500]);
    scatter(ec_wa(ok), ec_dmag(ok), 15, [0 0 0]+0.5, '+');
    hold on
    pS = scatter(ec_wa(~ok & Cinx), ec_dmag(~ok & Cinx), 400, ...
            ec_vmag(~ok & Cinx) + ec_dmag(~ok & Cinx), '.', gprops{:});
    ax = axis;

    style_wa_dmag_plot(...
        sprintf('%s: Earth Characterizations vs. WA and dMag, shaded by Planet Magnitude', full_name),...
        'Planet Apparent Visual Magnitude');
    legend({sprintf('Successful Chars (%d)', sum(ok)), ...
            sprintf('Failed %s Chars (%d)', full_name, sum(~ok & Cinx))}, ...
           'Location', 'northeast');
    % easy chars sub-legend
    text(txt_x, txt_y, { ...
        '"Easy" Characterizations within red box:', ...
        sprintf('  %d Successful', sum(easy & ok)), ...
        sprintf('  %d Failed (counting %s only)', sum(easy & Cinx & ~ok), full_name)}, ...
              'FontSize', 13);
    write_plots(sprintf('earth-char-wa-dmag-%s-vmag', abbr_name));

    % WA/dMag, shading = Adjusted Dmag
    if new_figs; figure(); end; clf;
    set(gcf, 'Position', [100 100 850 500]);
    % adjusted dmag - improvement possible
    ec_dmagPhi = ec_dmag + max(-5, 2.5*log10(min(0.7, ec_phi) / 0.7));
    % FIXME: dmag or dmagPhi for the OK's?
    scatter(ec_wa(ok), ec_dmag(ok), 15, [0 0 0]+0.5, '+');
    hold on
    % scatter(ec_wa(~ok & promo), ec_dmagPhi(~ok & promo), 40, ec_phi(~ok&promo),'s');
    pS = scatter(ec_wa(~ok & Cinx), ec_dmagPhi(~ok & Cinx), 400, ...
                 ec_phi(~ok&Cinx), '.', gprops{:});

    %set(pS, gprops{:});
    style_wa_dmag_plot(sprintf('%s: Earth Characterizations vs. WA and ADJUSTED dMag, shaded by Phi', full_name), ...
                       'Lambertian \Phi');
    legend({'Successful Chars', sprintf('Failed %s Chars', full_name)}, 'Location','northeast');
    write_plots(sprintf('earth-char-wa-dmag-%s-zzz-dmag', abbr_name));

end % for pnum


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Remaining plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% supplementary title
title_x = plot_make_title(t_info);

tprops = {'FontWeight', 'bold'};

%% histogram of Phi
if new_figs; figure(); end; clf;
set(gcf, 'Position', [100 100 850 500]);
colormap(brighten(jet, 0.5)); % stacked bars will be red/blue, like lines in plot
h_bins = [0:0.04:1];
h_promo_fail = histcounts(ec_phi(~ok & promo), h_bins);
h_deep_fail  = histcounts(ec_phi(~ok & deep),  h_bins);
h_promo_ok   = histcounts(ec_phi( ok & promo), h_bins);
h_deep_ok    = histcounts(ec_phi( ok & deep),  h_bins);
h_deep  = h_deep_fail + h_deep_ok;
h_promo = h_promo_fail + h_promo_ok;

% TODO: use "yyaxis right" instead of plotyy!

stakbar = @(x,y) bar(x, y, 'stacked');

h_bins1 = 0.5 * (h_bins(1:end-1) + h_bins(2:end));
[ax, h_left, h_right] = ...
    plotyy(h_bins1', [h_deep_fail; h_promo_fail]', ...
           h_bins1', [h_deep_fail./h_deep; h_promo_fail./h_promo], ...
           stakbar, @plot);
% bar(h_bins1', [h_deep; h_promo]', 'stacked');
%axis([0 1 0 Inf])
h_right(1).LineWidth = 2;
h_right(2).LineWidth = 2;

% style the plot
grid on
% title: prevent special interpretation of _
title1 = 'Failed Earth Characterizations vs. Lambertian Phi (Whole Ensemble)';
title({title_x, title1}, 'Interpreter', 'none');
xlabel('Lambertian Reflectance Phi', tprops{:});
ylabel(ax(1), 'Number of Failed Characterizations [count]', tprops{:});
ylabel(ax(2), 'Failure Probability', 'FontSize', 13, tprops{:}); % has to be given here
set(gca,'TitleFontSizeMultiplier', 1.1)
set(gca, 'FontSize', 13);
legend('Deep-dive Targets', 'Promoted Targets');
write_plots('earth-char-hist-phi');

%% histogram of Log Phi
if new_figs; figure(); end; clf;
set(gcf, 'Position', [100 100 850 500]);
colormap(brighten(jet, 0.5)); % stacked bars will be red/blue, like lines in plot
del = 0.4;
h_bins_log = [-6:del:0];
h_promo = histcounts(2.5*log10(ec_phi(~ok & promo)), [-Inf h_bins_log]);
h_deep  = histcounts(2.5*log10(ec_phi(~ok & deep)),  [-Inf h_bins_log]);
h_bins1 = 0.5 * (h_bins_log(1:end-1) + h_bins_log(2:end));
h_bins1 = [h_bins1(1)-del h_bins1]; % add leftmost bin
bar(h_bins1', [h_deep; h_promo]', 'stacked');

% style the plot
grid on;
% title: prevent special interpretation of _
title1 = 'Failed Earth Characterizations vs. Magnitude Difference (Whole Ensemble)';
title({title_x, title1}, 'Interpreter', 'none');
xlabel('Lambertian Reflectance Phi, as Magnitude (2.5 log_{10} \Phi)', tprops{:});
ylabel('Number of Failed Characterizations [count]', tprops{:});
legend('Deep-dive Targets', 'Promoted Targets');
set(gca,'TitleFontSizeMultiplier', 1.1)
set(gca, 'FontSize', 13);
write_plots('earth-char-hist-phi-log');

%keyboard

end
