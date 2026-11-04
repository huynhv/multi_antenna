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
% plot_as_matrix = true: one figure per attack, rows ROC / Pd vs threshold / Pf vs threshold.
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
plot_as_matrix = true;   % true: one figure per attack, 3 x 3 grid (rows: ROC, Pd vs threshold, Pf vs threshold;
                         % columns: tests). false: three figures (ROC, Pf vs threshold, Pd vs threshold), rows = attacks
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
tests = {"normcorr", R.roc_corr_grid, R.roc_corr_det ./ n_atk, R.roc_corr_fa ./ n_hon, ic, 'corr\_min'};
for param = ["alpha", "t0"]
    det = R.(char("roc_" + param + "_det"));   % k x spread x attack x channel x agent x sensor
    fa  = R.(char("roc_" + param + "_fa"));
    [~, ik] = min(abs(R.roc_k_grid - R.(char(param + "_delta_k"))));
    tests(end+1, :) = {param + " (Qn)", R.roc_k_grid, ...
        reshape(det(:, sp, :, :, :, :), [nk det_size]) ./ n_atk, ...
        reshape(fa(:, sp, :, :, :, :),  [nk det_size]) ./ n_hon, ik, 'k (Qn sigmas)'}; %#ok<SAGROW>
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

% Column headers: test name (* = used by the detector) and its current threshold
col_title = cell(n_tests, 1);
for s = 1:n_tests
    star = "";
    if in_detector(s)
        star = " *";
    end
    col_title{s} = {char(test_names(s) + star), sprintf('o: %s = %g', tests{s,6}, tests{s,2}(tests{s,5}))};
end
show_idx = 1:n_tests;
run_info = sprintf('agent SNR %g dB, S = %d, Pf source: %s, * = used by the detector', ...
    R.agent_db_values(agent_db_idx), R.sensor_vals(sensor_idx), pf_source);
fig_count = 0;

%% Matrix mode: one figure per attack, rows ROC / Pd vs threshold / Pf vs threshold, columns tests
% ("none" has no attacker sensors, so no figure; its Pf is in the tables below)
if plot_as_matrix
    row_names = ["ROC", "P_d vs threshold", "P_f vs threshold"];
    for a = atk_list
        fig_count = fig_count + 1;
        [tl, fig] = grid_figure(fig_count, 3, numel(show_idx), ...
            {sprintf('Attack: %s -- ROC, Pd and Pf vs threshold, one curve per channel SNR', R.attacks(a)), run_info});
        for r = 1:3
            for cc = 1:numel(show_idx)
                s = show_idx(cc);
                [~, grid_vals, Pd, Pf, i_now, xlab] = tests{s, :};
                ax = nexttile(tl);
                hold(ax, 'on');
                switch r
                    case 1
                        draw_roc(ax, Pd, Pf, i_now, a, agent_db_idx, sensor_idx, snr_vals, snr_colors, log_pf);
                        ylab = 'P_d';  xl = 'P_f';
                    case 2
                        draw_vs_thresh(ax, grid_vals, Pd, i_now, a, agent_db_idx, sensor_idx, snr_vals, snr_colors, false);
                        ylab = 'P_d';  xl = '';   % shares the x axis with row 3
                    case 3
                        draw_vs_thresh(ax, grid_vals, Pf, i_now, a, agent_db_idx, sensor_idx, snr_vals, snr_colors, log_pf);
                        ylab = 'P_f';  xl = xlab;
                end
                if r == 1
                    title(ax, col_title{s}, 'FontWeight', 'normal');
                end
                if cc == 1
                    ylabel(ax, {char(row_names(r)), ylab});
                elseif ~(r == 3 && log_pf)   % log-Pf tiles may not share y limits
                    ax.YTickLabel = [];
                end
                if isempty(xl)
                    ax.XTickLabel = [];
                else
                    xlabel(ax, xl);
                end
                if r == 1 && cc == 1
                    ax_leg = ax;
                end
            end
        end
        shared_legend(ax_leg, 'Channel SNR', show_legend);

        if save_figs
            home_dir = char(java.lang.System.getProperty('user.home'));
            folder_path = home_dir + "\Desktop\Byzantine AirComp Results\";

            save_folder = folder_path + experiment_name;
            if ~exist(save_folder, 'dir')
                mkdir(save_folder);
            end
            exportgraphics(fig, fullfile(save_folder, experiment_name + "_roc_matrix_" + R.attacks(a) + "_" + ...
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
fprintf('\nPf source: %s | agent SNR = %g dB | S = %d | detector uses: normcorr + %s\n', pf_source, ...
    R.agent_db_values(agent_db_idx), R.sensor_vals(sensor_idx), strjoin(test_names(in_detector & test_names ~= "normcorr"), ' + '));
for s = 1:n_tests
    print_table(sprintf('AUC -- %s (sweep %s)', tests{s,1}, tests{s,6}), auc(:, :, s), R.attacks, snr_vals);
end

% Operating point of each test at the thresholds used in the simulation (all attacks, incl. "none")
op = @(P, i) reshape(P(i, :, :, agent_db_idx, sensor_idx), n_atks, n_snr);
for s = 1:n_tests
    [name, grid_vals, Pd, Pf, i_now, xlab] = tests{s, :};
    print_table(sprintf('Pd -- %s at %s = %g', name, xlab, grid_vals(i_now)), op(Pd, i_now), R.attacks, snr_vals);
    print_table(sprintf('Pf -- %s at %s = %g', name, xlab, grid_vals(i_now)), op(Pf, i_now), R.attacks, snr_vals);
end


%% ---- Local functions ----
function [tl, fig] = grid_figure(fig_count, rows, cols, ttl)
% Figure sized to fit the grid on screen, with a tight tiled layout and an overall title.
    fig = figure("Name", ttl{1}, "Color", "w");
    scr = get(groot, 'ScreenSize');
    w = min(scr(3) - 80,  230*cols + 170);
    h = min(scr(4) - 140, 200*rows + 120);
    set(fig, 'Position', [40 + 25*fig_count, max(40, scr(4) - h - 90 - 25*fig_count), w, h]);
    tl = tiledlayout(fig, rows, cols, 'TileSpacing', 'tight', 'Padding', 'compact');
    title(tl, ttl{1}, 'FontWeight', 'bold');
    subtitle(tl, ttl{2});
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

function label_tile(ax, r, c, n_rows, row_name, col_header, xlab, ylab)
% Column headers on the top row, attack name + y label on the first column, x label on the
% bottom row; interior tick labels are hidden (all tiles in a row/column share limits).
    if r == 1
        title(ax, col_header, 'FontWeight', 'normal');
    end
    if c == 1
        ylabel(ax, {char(row_name), ylab});
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
