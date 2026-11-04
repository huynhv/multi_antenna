% PLOT_TIME_SHIFT  Time-shift sensitivity: MSE, Pd and Pf against the attacker's shift n*Delta.
% Run standalone, or after multi_antenna.m with experiment_list = "time_shift_sensitivity".
% The attacker's pulse is moved to t0_true + n*Delta, Delta = one sample (dt). x-axis ticks show n
% and the shift as a percentage of t0_true.
%   Figure 1: normalized MSE (alpha, t0) per estimator; dashed = proposed with no attack
%   Figure 2: Pd and Pf of each detection test at the thresholds used in the simulation
% A dashed vertical line marks the detection floor: the smallest shift with combined Pd >= pd_target.

if exist('experiment', 'var') && experiment == "time_shift_sensitivity"
    results_file = experiment + "_results.mat";
else
    results_file = "time_shift_sensitivity_results.mat";
end
R = load(results_file);   % struct, so the caller's workspace isn't overwritten
experiment_name = erase(string(results_file), "_results.mat");   % used in the saved file names

pd_target    = 0.9;   % detection floor = smallest shift with combined Pd >= pd_target
show_snr     = 5;    % channel SNR(s) (dB) to draw and tabulate, any of those in the run
save_figs    = true;  % save each figure as a 300 dpi PNG to ~\Desktop\Byzantine AirComp Results\<experiment>\
agent_db_idx = 1;
sensor_idx   = 1;
scheme_idx   = 1;

%% Sweep layout
n_vals   = R.ts_multiples(:).';                  % multiples of Delta
pct_vals = R.ts_shift_pct(:).';                  % shift as a percentage of t0_true
p_none   = find(R.pivot_types == "none", 1);     % no-attack reference
p_ts     = find(R.pivot_types == "time shift");  % one per n, in the order of ts_multiples
snr_vals = R.channel_db_values(scheme_idx, :, agent_db_idx);
n_snr    = numel(snr_vals);
snr_idx  = find(ismember(snr_vals, show_snr));   % channel SNRs drawn and tabulated
if isempty(snr_idx)
    error('None of show_snr (%s dB) is in these results (channel SNRs: %s dB).', mat2str(show_snr), mat2str(snr_vals));
end
norm_coeffs = [R.alpha_true^2, R.t0_true^2];     % normalized MSE = MSE / theta_true^2
param_lbl   = ["\alpha", "t_0"];
x_ticks  = n_vals;
x_labels = compose("%g\\newline(%.2g%%)", n_vals.', pct_vals.');
run_info = sprintf('\\Delta = dt = %.3g ms = %.3g%% of t_0 (%.4g s) | agent SNR %g dB, S = %d, detector: normcorr + %s', ...
    1e3*R.ts_delta, 100*R.ts_delta/R.t0_true, R.t0_true, R.agent_db_values(agent_db_idx), R.sensor_vals(sensor_idx), upper(R.detector_spread));

% MSE curves: {legend name, averaged MSE array, colour, line style}
C = orderedcolors("gem12");
methods = {
    'Proposed',    R.avg_dep_mse,             C(1,:), '-o'
    'Undefended',  R.avg_dep_mse_undefended,  C(4,:), '-s'
    'Oracle',      R.avg_dep_mse_oracle,      C(5,:), '-^'
    'Full honest', R.avg_dep_mse_full_honest, C(7,:), ':'
    };
conv_names = string(fieldnames(R.avg_dep_mse_conv)).';
conv_cols  = C([2 6 11 12], :);
for q = 1:numel(conv_names)
    methods(end+1, :) = {char(replace(conv_names(q), "_", " ")), R.avg_dep_mse_conv.(char(conv_names(q))), conv_cols(q,:), '--d'}; %#ok<SAGROW>
end
mse_of = @(A, prm, c, p) reshape(A(1, agent_db_idx, c, prm, scheme_idx, p, sensor_idx), 1, []) / norm_coeffs(prm);

% Detection rates at the simulation's thresholds: {legend name, Pd(p, c), Pf(p, c), colour}
[~, ic] = min(abs(R.roc_corr_grid - R.corr_min));
[~, ika] = min(abs(R.roc_k_grid - R.alpha_delta_k));
[~, ikt] = min(abs(R.roc_k_grid - R.t0_delta_k));
sp  = find(R.roc_spread_names == R.detector_spread, 1);
idx = @(p, c) sub2ind(size(R.roc_n_atk, 1:4), p, c*ones(size(p)), agent_db_idx*ones(size(p)), sensor_idx*ones(size(p)));
n_atk = @(p, c) R.roc_n_atk(idx(p, c));
n_hon = @(p, c) R.roc_n_hon(idx(p, c));
pick_k = @(A, ik, p, c) reshape(A(ik, sp, p, c, agent_db_idx, sensor_idx), 1, []);
tests = {
    'Combined detector', @(p,c) R.det_tp(idx(p,c)) ./ (R.det_tp(idx(p,c)) + R.det_fn(idx(p,c))), ...
                         @(p,c) R.det_fp(idx(p,c)) ./ n_hon(p,c), [0 0 0]
    sprintf('normcorr (corr_{min} = %g)', R.corr_min), ...
                         @(p,c) reshape(R.roc_corr_det(ic, p, c, agent_db_idx, sensor_idx), 1, []) ./ n_atk(p,c), ...
                         @(p,c) reshape(R.roc_corr_fa(ic, p, c, agent_db_idx, sensor_idx), 1, []) ./ n_hon(p,c), C(2,:)
    sprintf('t_0 test (k = %g)', R.t0_delta_k), ...
                         @(p,c) pick_k(R.roc_t0_det, ikt, p, c) ./ n_atk(p,c), @(p,c) pick_k(R.roc_t0_fa, ikt, p, c) ./ n_hon(p,c), C(1,:)
    sprintf('\\alpha test (k = %g)', R.alpha_delta_k), ...
                         @(p,c) pick_k(R.roc_alpha_det, ika, p, c) ./ n_atk(p,c), @(p,c) pick_k(R.roc_alpha_fa, ika, p, c) ./ n_hon(p,c), C(3,:)
    };

% Detection floor per channel SNR (combined detector)
floor_n = nan(1, n_snr);
for c = snr_idx
    pd_comb = tests{1,2}(p_ts, c);
    i_floor = find(pd_comb >= pd_target, 1);
    if ~isempty(i_floor)
        floor_n(c) = n_vals(i_floor);
    end
end

%% Figure 1: normalized MSE vs shift (rows: alpha, t0; columns: channel SNR)
[tl, fig1] = grid_figure(1, 2, numel(snr_idx), {'Time-shift sensitivity: normalized MSE vs shift n\Delta', run_info});
for prm = 1:2
    for c = snr_idx
        ax = nexttile(tl);
        hold(ax, 'on');
        for q = 1:size(methods, 1)
            plot(ax, n_vals, mse_of(methods{q,2}, prm, c, p_ts), methods{q,4}, 'Color', methods{q,3}, ...
                'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', methods{q,1});
        end
        yline(ax, mse_of(R.avg_dep_mse, prm, c, p_none), '--', 'Color', C(1,:), 'LineWidth', 1, ...
            'DisplayName', 'Proposed, no attack');
        mark_floor(ax, floor_n(c));
        ax.YScale = 'log';
        style_shift_axes(ax, x_ticks, x_labels);
        if prm == 1
            title(ax, sprintf('%g dB channel SNR', snr_vals(c)), 'FontWeight', 'normal');
        end
        if c == snr_idx(1)
            ylabel(ax, sprintf('NMSE(%s)', param_lbl(prm)));
        end
        if prm == 2
            xlabel(ax, 'n (shift as % of t_0)');
        end
        if prm == 1 && c == snr_idx(1)
            ax_leg = ax;
        end
    end
end
lgd = legend(ax_leg, 'Orientation', 'vertical');
lgd.Layout.Tile = 'east';

%% Figure 2: Pd and Pf vs shift (rows: Pd, Pf; columns: channel SNR)
[tl, fig2] = grid_figure(2, 2, numel(snr_idx), {sprintf('Time-shift sensitivity: detection vs shift n\\Delta (dashed vertical: Pd = %g floor)', pd_target), run_info});
for row = 1:2
    for c = snr_idx
        ax = nexttile(tl);
        hold(ax, 'on');
        for q = 1:size(tests, 1)
            rate = tests{q, 1 + row}(p_ts, c);
            plot(ax, n_vals, rate, '-o', 'Color', tests{q,4}, 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', tests{q,1});
        end
        if row == 1
            yline(ax, pd_target, ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
            ylim(ax, [-0.05 1.08]);
        else
            yline(ax, tests{1,3}(p_none, c), '--', 'Color', [0 0 0], 'LineWidth', 1, 'DisplayName', 'Combined Pf, no attack');
            pf_all = cellfun(@(f) max(f(p_ts, c)), tests(:,3));
            ylim(ax, [0, max(0.02, 1.15*max([pf_all; tests{1,3}(p_none, c)]))]);
        end
        mark_floor(ax, floor_n(c));
        style_shift_axes(ax, x_ticks, x_labels);
        if row == 1
            title(ax, sprintf('%g dB channel SNR', snr_vals(c)), 'FontWeight', 'normal');
        end
        if c == snr_idx(1)
            rate_lbl = ["P_d", "P_f"];
            ylabel(ax, rate_lbl(row));
        end
        if row == 2
            xlabel(ax, 'n (shift as % of t_0)');
        end
        if row == 2 && c == snr_idx(1)
            ax_leg = ax;   % the Pf row carries the extra "no attack" legend entry
        end
    end
end
lgd = legend(ax_leg, 'Orientation', 'vertical');
lgd.Layout.Tile = 'east';

if save_figs
    file_info = string(R.agent_db_values(agent_db_idx)) + "_db_Tp=" + string(R.Tp) + "_S=" + string(R.sensor_vals(sensor_idx));
    save_figure(fig1, experiment_name, "mse_" + file_info);
    save_figure(fig2, experiment_name, "detection_" + file_info);
end

%% Table: detection floor and worst-case harm per channel SNR
fprintf('\nTime-shift sensitivity (Delta = dt = %.4g ms = %.4g%% of t0_true)\n', 1e3*R.ts_delta, 100*R.ts_delta/R.t0_true);
fprintf('%9s | %-28s | %-34s | %-34s\n', 'SNR (dB)', sprintf('detection floor (Pd >= %g)', pd_target), ...
    'worst proposed/oracle NMSE, alpha', 'worst proposed/oracle NMSE, t0');
for c = snr_idx
    if isnan(floor_n(c))
        floor_str = 'not reached';
    else
        q = find(n_vals == floor_n(c), 1);
        floor_str = sprintf('n = %g (%.3g ms, %.3g%%)', floor_n(c), 1e3*floor_n(c)*R.ts_delta, pct_vals(q));
    end
    worst = strings(1, 2);
    for prm = 1:2
        ratio = mse_of(R.avg_dep_mse, prm, c, p_ts) ./ mse_of(R.avg_dep_mse_oracle, prm, c, p_ts);
        [rmax, q] = max(ratio);
        worst(prm) = sprintf('%.3g x at n = %g', rmax, n_vals(q));
    end
    fprintf('%9g | %-28s | %-34s | %-34s\n', snr_vals(c), floor_str, worst(1), worst(2));
end


%% ---- Local functions ----
function save_figure(fig, experiment_name, file_tag)
% Save a figure as a 300 dpi PNG to ~\Desktop\Byzantine AirComp Results\<experiment>\.
    home_dir = char(java.lang.System.getProperty('user.home'));
    folder_path = home_dir + "\Desktop\Byzantine AirComp Results\";

    save_folder = folder_path + experiment_name;
    if ~exist(save_folder, 'dir')
        mkdir(save_folder);
    end
    exportgraphics(fig, fullfile(save_folder, experiment_name + "_" + file_tag + ".png"), 'Resolution', 300);
end

function [tl, fig] = grid_figure(fig_count, rows, cols, ttl)
    fig = figure("Name", ttl{1}, "Color", "w");
    scr = get(groot, 'ScreenSize');
    w = min(scr(3) - 80,  300*cols + 220);
    h = min(scr(4) - 140, 280*rows + 140);
    set(fig, 'Position', [40 + 25*fig_count, max(40, scr(4) - h - 90 - 25*fig_count), w, h]);
    tl = tiledlayout(fig, rows, cols, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, ttl{1}, 'FontWeight', 'bold');
    subtitle(tl, ttl{2});
end

function style_shift_axes(ax, x_ticks, x_labels)
    ax.XScale = 'log';
    xlim(ax, [min(x_ticks)/1.3, max(x_ticks)*1.3]);
    ax.XTick = x_ticks;
    ax.XTickLabel = x_labels;
    grid(ax, 'on');
    box(ax, 'on');
end

function mark_floor(ax, n_floor)
    if ~isnan(n_floor)
        xline(ax, n_floor, '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 1, 'HandleVisibility', 'off');
    end
end
