% PLOT_TIME_SHIFT  Time-shift sensitivity: MSE, Pd and Pf against the attacker's shift n*Delta.
% Run standalone, or after multi_antenna.m with experiment_list = "time_shift_sensitivity".
% The attacker's pulse is moved to t0_true + n*Delta, Delta = one sample (dt). x-axis: n.
%   Figure 1: normalized MSE (alpha, t0) per estimator
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
show_snr     = 5;    % channel SNR(s) (dB) to draw and tabulate, any of those in the run (one row each)
leg_tile     = 2;    % Figure 1 tile that holds the legend: 1 = NMSE(alpha), 2 = NMSE(t0)
font_size    = 12;    % font size of all axes (legends: 2 pt smaller)
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
x_labels = compose("%g", n_vals.');
x_labels(2:2:end) = "";                          % label every other tick to avoid crowding
run_info = sprintf('\\Delta = dt = %.3g ms = %.3g%% of t_0 (%.4g s) | agent SNR %g dB, S = %d, detector: cosine similarity + %s', ...
    1e3*R.ts_delta, 100*R.ts_delta/R.t0_true, R.t0_true, R.agent_db_values(agent_db_idx), R.sensor_vals(sensor_idx), upper(R.detector_spread));

% MSE curves: {legend name, averaged MSE array, colour, line style}
C = orderedcolors("gem12");
methods = {
    'Proposed',    R.avg_dep_mse,             C(1,:), '-o'
    'Undefended',  R.avg_dep_mse_undefended,  C(4,:), '-s'
    'Oracle',      R.avg_dep_mse_oracle,      C(5,:), '-^'
    'Byzantine-Free', R.avg_dep_mse_full_honest, C(7,:), ':'
    };
conv_names = string(fieldnames(R.avg_dep_mse_conv)).';
conv_style = containers.Map({'median', 'trimmed_mean', 'skipped_mean', 'skipped_wmean'}, ...   % same as plot_attacks
    {{C(2,:), 'Median'}, {C(12,:), 'Trimmed mean'}, {C(6,:), 'Skipped mean'}, {C(11,:), 'Weighted skipped mean'}});
for q = 1:numel(conv_names)
    st = conv_style(char(conv_names(q)));   % {colour, legend name}
    methods(end+1, :) = {st{2}, R.avg_dep_mse_conv.(char(conv_names(q))), st{1}, '--d'}; %#ok<SAGROW>
end
mse_of = @(A, prm, c, p) reshape(A(1, agent_db_idx, c, prm, scheme_idx, p, sensor_idx), 1, []) / norm_coeffs(prm);

% Detection rates at the simulation's thresholds: {legend name, Pd(p, c), Pf(p, c), colour}
[~, ic] = min(abs(R.roc_corr_grid - R.corr_min));
[~, ika] = min(abs(R.roc_k_grid - R.alpha_delta_k));
[~, ikt] = min(abs(R.roc_k_grid - R.t0_delta_k));
sp  = find(R.roc_spread_names == R.detector_spread, 1);
idx = @(p, c) sub2ind(size(R.roc_n_atk, 1:4), p, c*ones(size(p)), agent_db_idx*ones(size(p)), sensor_idx*ones(size(p)));
% Rows throughout: indexing a vector-shaped array (e.g. a single channel SNR) keeps the array's
% orientation, so a column here would broadcast row ./ column into a matrix.
pick_lin = @(A, p, c) reshape(A(idx(p, c)), 1, []);
n_atk = @(p, c) pick_lin(R.roc_n_atk, p, c);
n_hon = @(p, c) pick_lin(R.roc_n_hon, p, c);
pick_k = @(A, ik, p, c) reshape(A(ik, sp, p, c, agent_db_idx, sensor_idx), 1, []);
tests = {
    'Combined detector', @(p,c) pick_lin(R.det_tp, p, c) ./ (pick_lin(R.det_tp, p, c) + pick_lin(R.det_fn, p, c)), ...
                         @(p,c) pick_lin(R.det_fp, p, c) ./ n_hon(p,c), [0 0 0]
    'Cosine similarity', ...
                         @(p,c) reshape(R.roc_corr_det(ic, p, c, agent_db_idx, sensor_idx), 1, []) ./ n_atk(p,c), ...
                         @(p,c) reshape(R.roc_corr_fa(ic, p, c, agent_db_idx, sensor_idx), 1, []) ./ n_hon(p,c), C(2,:)
    't_0 test', ...
                         @(p,c) pick_k(R.roc_t0_det, ikt, p, c) ./ n_atk(p,c), @(p,c) pick_k(R.roc_t0_fa, ikt, p, c) ./ n_hon(p,c), C(1,:)
    '\alpha test', ...
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

%% Figure 1: normalized MSE vs shift (columns: alpha, t0; rows: channel SNR)
[tl, fig1] = grid_figure(1, numel(snr_idx), 2, {'Time-shift sensitivity: normalized MSE vs shift n\Delta', run_info});
for c = snr_idx
    for prm = 1:2
        ax = nexttile(tl);
        hold(ax, 'on');
        for q = 1:size(methods, 1)
            plot(ax, n_vals, mse_of(methods{q,2}, prm, c, p_ts), methods{q,4}, 'Color', methods{q,3}, ...
                'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', methods{q,1});
        end
        mark_floor(ax, floor_n(c));
        ax.YScale = 'log';
        style_shift_axes(ax, x_ticks, x_labels, font_size);
        if c == snr_idx(1)
            title(ax, sprintf('NMSE(%s)', param_lbl(prm)), 'FontWeight', 'normal');
        end
        if prm == 1
            ylabel(ax, 'Normalized MSE');
        end
        if c == snr_idx(end)
            xlabel(ax, 'n  (attacker delay t_0 + n\Delta{}t)');
        end
        if prm == leg_tile && c == snr_idx(1)
            legend(ax, 'Location', 'best', 'FontSize', font_size - 2);
        end
    end
end

%% Figure 2: Pd and Pf vs shift (columns: Pd, Pf; rows: channel SNR)
[tl, fig2] = grid_figure(2, numel(snr_idx), 2, {sprintf('Time-shift sensitivity: detection vs shift n\\Delta (dashed vertical: Pd = %g floor)', pd_target), run_info});
for c = snr_idx
    for row = 1:2   % 1: Pd, 2: Pf (one column each)
        ax = nexttile(tl);
        hold(ax, 'on');
        for q = 1:size(tests, 1)
            rate = tests{q, 1 + row}(p_ts, c);
            plot(ax, n_vals, rate, '-o', 'Color', tests{q,4}, 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', tests{q,1});
        end
        if row == 1
            yline(ax, pd_target, ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
            ylim(ax, [-0.05 1.08]);
            set_ticks(ax, 'y', 0:0.25:1);
        else
            pf_all = cellfun(@(f) max(f(p_ts, c)), tests(:,3));
            ylim(ax, [0, max(0.02, 1.15*max(pf_all))]);
        end
        mark_floor(ax, floor_n(c));
        style_shift_axes(ax, x_ticks, x_labels, font_size);
        rate_lbl = ["P_d", "P_f"];
        if c == snr_idx(1)
            title(ax, rate_lbl(row), 'FontWeight', 'normal');
        end
        if row == 1
            ylabel(ax, 'Rate');
        end
        if c == snr_idx(end)
            xlabel(ax, 'n  (attacker delay t_0 + n\Delta{}t)');
        end
        if row == 2 && c == snr_idx(1)   % legend in the Pf tile
            legend(ax, 'Location', 'best', 'FontSize', font_size - 2);
        end
    end
end

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
    w = min(scr(3) - 80,  380*cols + 100);
    h = min(scr(4) - 140, 280*rows + 140);
    set(fig, 'Position', [40 + 25*fig_count, max(40, scr(4) - h - 90 - 25*fig_count), w, h]);
    tl = tiledlayout(fig, rows, cols, 'TileSpacing', 'compact', 'Padding', 'compact');
end

function style_shift_axes(ax, x_ticks, x_labels, font_size)
    ax.FontSize = font_size;
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
