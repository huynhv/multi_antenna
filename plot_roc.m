% PLOT_ROC  Per-sensor Pd vs Pf for each detection test, with each test's own threshold swept.
% Run standalone, or after multi_antenna.m (which defines `experiment`).
%   Pd = fraction of (attacker sensor, trial, deployment) cases flagged
%   Pf = fraction of (honest sensor, trial, deployment) cases flagged
%   normcorr : sweep corr_min (flagged if norm_corr < corr_min)
%   alpha/t0 : sweep k of the robust outlier test, flagged if |x_i - median| / Qn > k
% Three tests (columns): normcorr, alpha (Qn), t0 (Qn).
% plot_as_matrix = false: three figures, each a grid with one row per attack and one column per test:
%   1. ROC (Pd vs Pf), one curve per channel SNR, 'o' at the threshold used in the simulation
%   2. Pf vs threshold        3. Pd vs threshold
% plot_as_matrix = true: one figure per attack, 2 x 2. Top: the combined detector's Pd / Pf / per-trial
%   rates vs channel SNR, then the ROC. Bottom: Pd vs threshold, Pf vs threshold. ROC and threshold
%   tiles are at one channel SNR with the three tests overlaid (threshold tiles: corr_min on the
%   bottom x axis, k on the top x axis).
if exist('experiment', 'var')
    results_file = experiment + "_results.mat";
else
    results_file = "multi_antenna_base_results.mat";
end
experiment_name = erase(string(results_file), "_results.mat");   % used in the saved file names
R = load(results_file);   % struct, so the caller's workspace isn't overwritten

pf_source    = "same";   % "same": Pf from honest sensors in the same attack scenario
                         % "none": Pf from the no-attack scenario (Neyman-Pearson calibration; needs "none" in attacks)
log_pf       = false;    % true: log-scale Pf axis -- x on the ROC plots, y on the Pf-vs-threshold plots
                         % (points with Pf = 0 are then not drawn)
show_legend  = true;
plot_as_matrix = true;   % true: one figure per attack, 2 x 2 tiles: combined detector vs channel SNR, ROC (top);
                         % Pd vs threshold, Pf vs threshold (bottom); tests overlaid, at matrix_snr. false: three figures (ROC, Pf vs threshold, Pd vs
                         % threshold), rows = attacks, columns = tests, one curve per channel SNR
matrix_snr   = 5;        % channel SNR (dB) of the matrix figures
matrix_font  = 12;       % font size of the matrix figures (legends: 2 pt smaller)
save_figs    = true;     % matrix mode only: save each figure as a 300 dpi PNG to
                         % ~\Desktop\Byzantine AirComp Results\<experiment>\
agent_db_idx = 1;
sensor_idx   = 1;
scheme_idx   = 1;

%% Rates: counts / number of cases, one entry per test
% {name, threshold grid, Pd, Pf, index of the current threshold, threshold label}
% Pd / Pf arrays: threshold x attack x channel x agent x sensor
det_size = size(R.roc_n_atk, 1:4);
n_atk = reshape(R.roc_n_atk, [1 det_size]);
n_hon = reshape(R.roc_n_hon, [1 det_size]);
n_atks = numel(R.attacks);
nk = numel(R.roc_k_grid);
sp = find(R.roc_spread_names == "qn", 1);   % alpha / t0 tests use the Qn spread

[~, ic] = min(abs(R.roc_corr_grid - R.corr_min));
tests = {"Cosine similarity", R.roc_corr_grid, R.roc_corr_det ./ n_atk, R.roc_corr_fa ./ n_hon, ic, '\eta_{cs}'};
for param = ["t0", "alpha"]   % same test order as the other plot scripts: cosine similarity, t0, alpha
    det = R.(char("roc_" + param + "_det"));   % k x spread x attack x channel x agent x sensor
    fa  = R.(char("roc_" + param + "_fa"));
    [~, ik] = min(abs(R.roc_k_grid - R.(char(param + "_delta_k"))));
    tests(end+1, :) = {replace(param, ["alpha", "t0"], ["\alpha", "t_0"]) + " test", R.roc_k_grid, ...   % TeX names
        reshape(det(:, sp, :, :, :, :), [nk det_size]) ./ n_atk, ...
        reshape(fa(:, sp, :, :, :, :),  [nk det_size]) ./ n_hon, ik, '\kappa_{x}'}; %#ok<SAGROW>
end
n_tests = size(tests, 1);
test_names = string(tests(:, 1));
in_detector = [true; R.detector_spread == "qn"; R.detector_spread == "qn"];   % tests used for nulling

if pf_source == "none"
    none_idx = find(R.attacks == "none", 1);
    if isempty(none_idx)
        error('pf_source = "none" needs "none" in the attacks list of the run.');
    end
    for s = 1:n_tests
        tests{s, 4} = repmat(tests{s, 4}(:, none_idx, :, :, :), [1 n_atks 1 1 1]);
    end
end

snr_vals = R.channel_db_values(scheme_idx, :, agent_db_idx);
n_snr = numel(snr_vals);
snr_colors = parula(n_snr + 1);   % drop the last (yellow) entry for visibility
atk_list = find(R.roc_n_atk(:, 1, agent_db_idx, sensor_idx) > 0).';   % attacks with attacker sensors ("none" has no Pd)

% AUC of every test
auc = nan(n_atks, n_snr, n_tests);
for s = 1:n_tests
    for a = atk_list
        for c = 1:n_snr
            auc(a, c, s) = roc_auc(tests{s,4}(:, a, c, agent_db_idx, sensor_idx), tests{s,3}(:, a, c, agent_db_idx, sensor_idx));
        end
    end
end

% Column headers: test name and its current threshold
col_title = cell(n_tests, 1);
for s = 1:n_tests
    col_title{s} = {char(test_names(s)), sprintf('o: %s = %g', tests{s,6}, tests{s,2}(tests{s,5}))};
end
show_idx = 1:n_tests;
run_info = sprintf('agent SNR %g dB, S = %d, Pf source: %s', ...
    R.agent_db_values(agent_db_idx), R.sensor_vals(sensor_idx), pf_source);
fig_count = 0;

%% Matrix mode: one figure per attack, 2 x 2 tiles. Top row: the combined detector (OR of the tests at
% the simulation's thresholds) vs channel SNR, then the ROC at matrix_snr. Bottom row: Pd and Pf vs
% threshold at matrix_snr. ROC and threshold tiles overlay the three tests; threshold tiles have two
% x axes: corr_min (normcorr, bottom) and k (alpha / t0, top). ("none" has no attacker sensors, so no figure.)
if plot_as_matrix
    c = find(snr_vals == matrix_snr, 1);
    if isempty(c)
        error('matrix_snr = %g dB is not in these results (channel SNRs: %s dB).', matrix_snr, mat2str(snr_vals));
    end
    test_colors = orderedcolors("gem12");
    test_colors = test_colors([2 1 3], :);   % cosine similarity orange, t0 blue, alpha yellow (as in the other plot scripts)
    test_styles = ["-", "--", "-."];
    k_tests = 2:3;   % tests swept over k (top axis); test 1 (normcorr) is swept over corr_min (bottom axis)
    % Tick positions; set_ticks labels every other one (first, third, ...) to avoid crowding at matrix_font
    prob_ticks = 0:0.25:1;                                             % Pd / Pf / rate axes
    corr_ticks = 0:0.25:1;                                             % corr_min
    k_ticks    = 0:1:max(R.roc_k_grid(isfinite(R.roc_k_grid)));        % k (Qn sigmas)
    for a = atk_list
        fig_count = fig_count + 1;
        [tl, fig] = grid_figure(fig_count, 2, 2, ...
            {sprintf('Attack: %s -- %g dB channel SNR', R.attacks(a), matrix_snr), run_info}, false);   % no overall title: more room for the curves
        fig.Position(3:4) = [950 820];   % room for the second (top) x axes

        % Tile 1 -- combined detector (normcorr OR alpha OR t0 at the simulation's thresholds) vs channel SNR
        pf_a = a;   % attack scenario the false-alarm rates come from
        if pf_source == "none"
            pf_a = none_idx;
        end
        cnt = @(fld, aa) reshape(R.(fld)(aa, :, agent_db_idx, sensor_idx), 1, []);
        pd_comb  = cnt('det_tp', a) ./ (cnt('det_tp', a) + cnt('det_fn', a));
        pf_comb  = cnt('det_fp', pf_a) ./ cnt('roc_n_hon', pf_a);
        all_cght = cnt('det_all_caught', a) ./ cnt('det_ntrials', a);
        any_fa   = cnt('det_any_fa', pf_a) ./ cnt('det_ntrials', pf_a);
        ax = nexttile(tl);
        hold(ax, 'on');
        det_color = [0 0 0];
        fa_color  = [0.85 0.2 0.2];
        plot(ax, snr_vals, pd_comb,  '-o',  'Color', det_color, 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', 'P_d (per sensor)');
        plot(ax, snr_vals, all_cght, '--s', 'Color', det_color, 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', 'P(all attackers flagged)');
        plot(ax, snr_vals, pf_comb,  '-o',  'Color', fa_color,  'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', 'P_f (per sensor)');
        plot(ax, snr_vals, any_fa,   '--s', 'Color', fa_color,  'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', 'P(\geq 1 honest flagged)');
        plot(ax, matrix_snr*[1 1], [0 1], ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');   % SNR of the other tiles
        style_axes(ax, [min(snr_vals) max(snr_vals)], false, false);
        ax.FontSize = matrix_font;
        set_ticks(ax, 'x', snr_vals);
        set_ticks(ax, 'y', prob_ticks);
        title(ax, 'Combined detector vs Channel SNR', 'FontWeight', 'normal');
        xlabel(ax, 'Channel SNR (dB)');
        ylabel(ax, 'Rate');
        if show_legend
            legend(ax, 'Location', 'best', 'FontSize', matrix_font - 2);
        end

        % Tile 2 -- ROC: all tests on one axes, 'o' at the threshold used in the simulation
        ax = nexttile(tl);
        hold(ax, 'on');
        for s = 1:n_tests
            [~, ~, Pd, Pf, i_now, ~] = tests{s, :};
            pd = Pd(:, a, c, agent_db_idx, sensor_idx);
            pf = Pf(:, a, c, agent_db_idx, sensor_idx);
            plot(ax, pf, pd, test_styles(s), 'Color', test_colors(s,:), 'LineWidth', 1.5, 'DisplayName', col_title{s}{1});
            plot(ax, pf(i_now), pd(i_now), 'o', 'Color', test_colors(s,:), 'MarkerFaceColor', test_colors(s,:), 'HandleVisibility', 'off');
        end
        plot(ax, [0 1], [0 1], ':', 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off');   % chance line
        style_axes(ax, [0 1], log_pf, false);
        ax.FontSize = matrix_font;
        if ~log_pf
            set_ticks(ax, 'x', prob_ticks);
        end
        set_ticks(ax, 'y', prob_ticks);
        title(ax, 'ROC', 'FontWeight', 'normal');
        xlabel(ax, 'P_f');
        ylabel(ax, 'P_d');
        if show_legend
            legend(ax, 'Location', 'southeast', 'FontSize', matrix_font - 2);
        end

        % Tiles 3 and 4 -- Pd and Pf vs threshold, side by side
        for metric = ["P_d", "P_f"]
            log_y = log_pf && metric == "P_f";
            tile = 3 + (metric == "P_f");
            % bottom axes: normcorr vs corr_min
            [~, grid_vals, Pd, Pf, i_now, xlab] = tests{1, :};
            P = Pd;
            if metric == "P_f"
                P = Pf;
            end
            ax = nexttile(tl);
            hold(ax, 'on');
            fin = isfinite(grid_vals);   % drop the +/-Inf closure points
            plot(ax, grid_vals(fin), reshape(P(fin, a, c, agent_db_idx, sensor_idx), [], 1), test_styles(1), ...
                'Color', test_colors(1,:), 'LineWidth', 1.5, 'DisplayName', col_title{1}{1});
            plot(ax, grid_vals(i_now)*[1 1], [0 1], ':', 'Color', test_colors(1,:), 'HandleVisibility', 'off');
            style_axes(ax, [min(grid_vals(fin)) max(grid_vals(fin))], false, log_y);
            ax.FontSize = matrix_font;
            set_ticks(ax, 'x', corr_ticks);
            if ~log_y
                set_ticks(ax, 'y', prob_ticks);
            end
            ax.XColor = test_colors(1,:);
            ax.Box = 'off';
            xlabel(ax, '$\eta_{\mathrm{cs}}$', 'Interpreter', 'latex');   % normcorr threshold
            ylabel(ax, char(metric));

            % top axes, same tile: alpha / t0 tests vs k
            ax2 = axes(tl);
            ax2.Layout.Tile = tile;
            hold(ax2, 'on');
            for s = k_tests
                [~, grid_vals, Pd, Pf, i_now, xlab] = tests{s, :};
                P = Pd;
                if metric == "P_f"
                    P = Pf;
                end
                fin = isfinite(grid_vals);
                plot(ax2, grid_vals(fin), reshape(P(fin, a, c, agent_db_idx, sensor_idx), [], 1), test_styles(s), ...
                    'Color', test_colors(s,:), 'LineWidth', 1.5, 'DisplayName', col_title{s}{1});
                plot(ax2, grid_vals(i_now)*[1 1], [0 1], ':', 'Color', test_colors(s,:), 'HandleVisibility', 'off');
            end
            style_axes(ax2, [min(grid_vals(fin)) max(grid_vals(fin))], false, log_y);
            ax2.FontSize = matrix_font;
            set_ticks(ax2, 'x', k_ticks);
            if ~log_y
                ax2.YTick = prob_ticks;   % same y ticks as the bottom axes (labels hidden below)
            end
            grid(ax2, 'off');   % the bottom axes carries the grid
            ax2.XAxisLocation = 'top';
            ax2.YAxisLocation = 'right';
            ax2.YTickLabel = [];
            ax2.Color = 'none';
            ax2.Box = 'off';
            xlabel(ax2, '$\kappa_{x}$ ($\alpha$, $t_0$)', 'Interpreter', 'latex');   % alpha / t0 test threshold
            title(ax2, char(metric + " vs threshold"), 'FontWeight', 'normal');   % on the top axes, above its labels
            if log_y   % both axes share the y scale (linear limits are already equal, from style_axes)
                p_all = cellfun(@(T) reshape(T{4}(isfinite(T{2}), a, c, agent_db_idx, sensor_idx), [], 1), ...
                    num2cell(tests, 2), 'UniformOutput', false);
                p_all = vertcat(p_all{:});
                y_lim = [min(p_all(p_all > 0)) / 2, 1.1];
                ylim(ax, y_lim);
                ylim(ax2, y_lim);
            end
        end

        if save_figs
            home_dir = char(java.lang.System.getProperty('user.home'));
            folder_path = home_dir + "\Desktop\Byzantine AirComp Results\";

            save_folder = folder_path + experiment_name;
            if ~exist(save_folder, 'dir')
                mkdir(save_folder);
            end
            exportgraphics(fig, fullfile(save_folder, experiment_name + "_roc_matrix_" + R.attacks(a) + "_channel_" + ...
                string(matrix_snr) + "_db_agent_" + ...
                string(R.agent_db_values(agent_db_idx)) + "_db_Tp=" + string(R.Tp) + "_S=" + string(R.sensor_vals(sensor_idx)) + ...
                "_pf_" + pf_source + ".png"), 'Resolution', 300);
        end
    end
end

%% Figure 1: ROC grid (rows: attacks, columns: tests), one curve per channel SNR
if ~plot_as_matrix && ~isempty(show_idx) && ~isempty(atk_list)
    fig_count = fig_count + 1;
    tl = grid_figure(fig_count, numel(atk_list), numel(show_idx), ...
        {'ROC per test (Pd vs Pf), one curve per channel SNR', run_info});
    for r = 1:numel(atk_list)
        a = atk_list(r);
        for cc = 1:numel(show_idx)
            s = show_idx(cc);
            [~, ~, Pd, Pf, i_now, ~] = tests{s, :};
            ax = nexttile(tl);
            hold(ax, 'on');
            draw_roc(ax, Pd, Pf, i_now, a, agent_db_idx, sensor_idx, snr_vals, snr_colors, log_pf);
            label_tile(ax, r, cc, numel(atk_list), R.attacks(a), col_title{s}, 'P_f', 'P_d');
            if r == 1 && cc == 1
                ax_leg = ax;
            end
        end
    end
    shared_legend(ax_leg, 'Channel SNR', show_legend);
end

%% Figures 2 and 3: Pf and Pd vs each test's threshold (rows: attacks, columns: tests)
% Pf rows include "none"; Pd rows only attacks with attacker sensors.
for metric = ["P_f", "P_d"]
    if metric == "P_f"
        show_atks = 1:n_atks;
        fig_ttl = {'Pf vs threshold (dashed: current value), one curve per channel SNR', run_info};
    else
        show_atks = atk_list;
        fig_ttl = {'Pd vs threshold (dashed: current value), one curve per channel SNR', run_info};
    end
    if plot_as_matrix || isempty(show_idx) || isempty(show_atks)
        continue
    end
    fig_count = fig_count + 1;
    tl = grid_figure(fig_count, numel(show_atks), numel(show_idx), fig_ttl);
    for r = 1:numel(show_atks)
        a = show_atks(r);
        for cc = 1:numel(show_idx)
            s = show_idx(cc);
            [~, grid_vals, Pd, Pf, i_now, xlab] = tests{s, :};
            if metric == "P_f"
                P = Pf;
            else
                P = Pd;
            end
            ax = nexttile(tl);
            hold(ax, 'on');
            draw_vs_thresh(ax, grid_vals, P, i_now, a, agent_db_idx, sensor_idx, snr_vals, snr_colors, log_pf && metric == "P_f");
            label_tile(ax, r, cc, numel(show_atks), R.attacks(a), col_title{s}, xlab, char(metric));
            if r == 1 && cc == 1
                ax_leg = ax;
            end
        end
    end
    shared_legend(ax_leg, 'Channel SNR', show_legend);
end

%% Tables
plain = @(s) regexprep(s, '[\\{}]', '');   % TeX name -> console text (\eta_{cs} -> eta_cs, \alpha -> alpha)
fprintf('\nPf source: %s | agent SNR = %g dB | S = %d | detector uses: %s + %s\n', pf_source, ...
    R.agent_db_values(agent_db_idx), R.sensor_vals(sensor_idx), test_names(1), plain(strjoin(test_names(in_detector & (1:n_tests).' > 1), ' + ')));
for s = 1:n_tests
    print_table(sprintf('AUC -- %s (sweep %s)', plain(tests{s,1}), plain(tests{s,6})), auc(:, :, s), R.attacks, snr_vals);
end

% Operating point of each test at the thresholds used in the simulation (all attacks, incl. "none")
op = @(P, i) reshape(P(i, :, :, agent_db_idx, sensor_idx), n_atks, n_snr);
for s = 1:n_tests
    [name, grid_vals, Pd, Pf, i_now, xlab] = tests{s, :};
    print_table(sprintf('Pd -- %s at %s = %g', plain(name), plain(xlab), grid_vals(i_now)), op(Pd, i_now), R.attacks, snr_vals);
    print_table(sprintf('Pf -- %s at %s = %g', plain(name), plain(xlab), grid_vals(i_now)), op(Pf, i_now), R.attacks, snr_vals);
end


%% ---- Local functions ----
function [tl, fig] = grid_figure(fig_count, rows, cols, ttl, show_title)
% Figure sized to fit the grid on screen, with a tight tiled layout and (unless show_title is false)
% an overall title and subtitle. ttl{1} also names the figure window.
    if nargin < 5
        show_title = true;
    end
    fig = figure("Name", ttl{1}, "Color", "w");
    scr = get(groot, 'ScreenSize');
    w = min(scr(3) - 80,  230*cols + 170);
    h = min(scr(4) - 140, 200*rows + 120);
    set(fig, 'Position', [40 + 25*fig_count, max(40, scr(4) - h - 90 - 25*fig_count), w, h]);
    tl = tiledlayout(fig, rows, cols, 'TileSpacing', 'tight', 'Padding', 'compact');
    if show_title
        title(tl, ttl{1}, 'FontWeight', 'bold');
        subtitle(tl, ttl{2});
    end
end

function draw_roc(ax, Pd, Pf, i_now, a, agent_db_idx, sensor_idx, snr_vals, snr_colors, log_pf)
% ROC tile: Pd vs Pf, one curve per channel SNR, 'o' at the current threshold, chance line.
    for c = 1:numel(snr_vals)
        pd = Pd(:, a, c, agent_db_idx, sensor_idx);
        pf = Pf(:, a, c, agent_db_idx, sensor_idx);
        plot(ax, pf, pd, '-', 'Color', snr_colors(c,:), 'LineWidth', 1.5, 'DisplayName', sprintf('%g dB', snr_vals(c)));
        plot(ax, pf(i_now), pd(i_now), 'o', 'Color', snr_colors(c,:), 'MarkerFaceColor', snr_colors(c,:), 'HandleVisibility', 'off');
    end
    plot(ax, [0 1], [0 1], ':', 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off');   % chance line
    style_axes(ax, [0 1], log_pf, false);
end

function draw_vs_thresh(ax, grid_vals, P, i_now, a, agent_db_idx, sensor_idx, snr_vals, snr_colors, log_y)
% Rate vs threshold tile, one curve per channel SNR, dashed line at the current threshold.
    fin = isfinite(grid_vals);   % drop the +/-Inf closure points
    xv = grid_vals(fin);
    for c = 1:numel(snr_vals)
        plot(ax, xv, reshape(P(fin, a, c, agent_db_idx, sensor_idx), [], 1), '-', 'Color', snr_colors(c,:), ...
            'LineWidth', 1.5, 'DisplayName', sprintf('%g dB', snr_vals(c)));
    end
    plot(ax, grid_vals(i_now)*[1 1], [0 1], '--', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
    style_axes(ax, [min(xv) max(xv)], false, log_y);
end

function set_ticks(ax, which, ticks)
% Ticks at every value of `ticks`, labels on every other one (the first, third, ...).
    lbl = cellstr(compose("%g", ticks(:)));
    lbl(2:2:end) = {''};
    if which == 'x'
        ax.XTick = ticks;
        ax.XTickLabel = lbl;
    else
        ax.YTick = ticks;
        ax.YTickLabel = lbl;
    end
end

function style_axes(ax, x_limits, log_x, log_y)
    grid(ax, 'on');
    box(ax, 'on');
    if log_x
        ax.XScale = 'log';
    elseif isequal(x_limits, [0 1])
        xlim(ax, [-0.05 1.05]);   % probability axis (ROC Pf): same margin as y
    else
        xlim(ax, x_limits);
    end
    if log_y
        ax.YScale = 'log';
    else
        ylim(ax, [-0.05 1.08]);   % margin so curves sitting at exactly 0 or 1 stay visible
    end
end

function lbl = attack_label(atk)
% Attack name as shown in the figures (others: title case, "time shift" -> "Time Shift")
    switch string(atk)
        case "noise"
            lbl = "Wideband Jamming";
        case "sinusoid"
            lbl = "Narrowband Jamming";
        otherwise
            lbl = string(regexprep(char(atk), '(^|\s)(\w)', '$1${upper($2)}'));
    end
end

function label_tile(ax, r, c, n_rows, row_name, col_header, xlab, ylab)
% Column headers on the top row, attack name + y label on the first column, x label on the
% bottom row; interior tick labels are hidden (all tiles in a row/column share limits).
    if r == 1
        title(ax, col_header, 'FontWeight', 'normal');
    end
    if c == 1
        ylabel(ax, {char(attack_label(row_name)), ylab});
    else
        ax.YTickLabel = [];
    end
    if r == n_rows
        xlabel(ax, xlab);
    else
        ax.XTickLabel = [];
    end
end

function shared_legend(ax, ttl, show)
% One legend for the whole grid, to the right of the tiles.
    if ~show
        return
    end
    lgd = legend(ax, 'Orientation', 'vertical');
    lgd.Layout.Tile = 'east';
    lgd.Title.String = ttl;
end

function a = roc_auc(pf, pd)
% Trapezoidal area under (pf, pd), closed with (0,0) and (1,1). Straight segments between
% grid points correspond to randomizing between neighbouring thresholds.
    pts = sortrows([0 0; pf(:) pd(:); 1 1]);
    pts = pts(all(isfinite(pts), 2), :);
    a = trapz(pts(:,1), pts(:,2));
end

function print_table(ttl, X, attacks, snr_vals)
    fprintf('\n-- %s --\n', ttl);
    disp(array2table(X, 'RowNames', cellstr(attacks), 'VariableNames', compose('%g dB', snr_vals)))
end
