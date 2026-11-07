% PLOT_CSI  Imperfect-CSI sensitivity: estimation and detection against the server's CSI NMSE.
% Run standalone, or after multi_antenna.m with experiment_list = "csi_sensitivity". Loads every
% csi_sensitivity_*_results.mat in the current folder (one file per CSI level); perfect CSI
% (NMSE = -Inf) is drawn at the left end of the x-axis and labelled "perfect".
%   Figure 1: normalized MSE (columns t0, alpha) vs CSI NMSE for the estimators in show_methods, one row
%             per channel SNR in show_snr.
%             A dashed vertical line marks the worst CSI at which the proposed t0 NMSE is still within
%             oracle_tol of the oracle.
%   Figure 2: detection vs CSI NMSE (columns Pd, Pf), one row per channel SNR in show_snr;
%             Pd and Pf per test at the thresholds used in the simulation.

files = dir("csi_sensitivity_*_results.mat");
if isempty(files)
    error('No csi_sensitivity_*_results.mat files in %s. Run multi_antenna.m with experiment_list = "csi_sensitivity".', pwd);
end
L = arrayfun(@(f) load(fullfile(f.folder, f.name)), files, 'UniformOutput', false);
csi_db = cellfun(@(x) x.csi_nmse_db, L);
[csi_db, order] = sort(csi_db);
L = L(order);

attack       = "noise";   % attack to plot (must be in csi_attacks of the run)
oracle_tol   = 0.10;      % "near-oracle": proposed t0 NMSE <= (1 + oracle_tol) x oracle
log_pf       = false;     % true: log-scale Pf axis (points with Pf = 0 are then not drawn)
show_snr     = 5;         % channel SNR(s) (dB) to draw and tabulate, any of those in the run (one row each)
leg_tile     = 1;         % Figure 1 tile that holds the legend: 1 = NMSE(t0), 2 = NMSE(alpha)
font_size    = 12;        % font size of all axes (legends: 2 pt smaller)
save_figs    = true;      % save each figure as a 300 dpi PNG to ~\Desktop\Byzantine AirComp Results\csi_sensitivity\
show_methods = ["Proposed", "Oracle", "median"];   % estimators drawn in Figure 1, any of: "Proposed",
                          % "Undefended", "Oracle", "Full honest", and the conventional baselines of the run
                          % ("median", "trimmed mean", "skipped mean", "skipped wmean")
agent_db_idx = 1;
sensor_idx   = 1;
scheme_idx   = 1;

%% Sweep layout
R1 = L{1};
a_idx = find(R1.attacks == attack, 1);
if isempty(a_idx)
    error('Attack "%s" is not in these results (attacks: %s).', attack, strjoin(R1.attacks, ', '));
end
snr_vals = R1.channel_db_values(scheme_idx, :, agent_db_idx);
n_snr    = numel(snr_vals);
snr_idx  = find(ismember(snr_vals, show_snr));   % channel SNRs drawn and tabulated
if isempty(snr_idx)
    error('None of show_snr (%s dB) is in these results (channel SNRs: %s dB).', mat2str(show_snr), mat2str(snr_vals));
end
n_lev   = numel(L);
norm_coeffs = [R1.alpha_true^2, R1.t0_true^2];   % normalized MSE = MSE / theta_true^2
param_lbl   = ["\alpha", "t_0"];

% x positions: finite NMSE levels as they are; perfect CSI one step to the left of the best level
fin = isfinite(csi_db);
step = 5;
if nnz(fin) > 1
    step = min(diff(csi_db(fin)));
end
x_pos = csi_db;
x_pos(~fin) = min(csi_db(fin)) - step;
x_labels = compose("%g", csi_db);
x_labels(~fin) = "perfect";
run_info = sprintf('attack: %s, compromised %s of S = %d | agent SNR %g dB | %d deployments x %d trials | detector: normcorr + %s', ...
    attack, mat2str(R1.pivot_atk{a_idx}), R1.sensor_vals(sensor_idx), R1.agent_db_values(agent_db_idx), ...
    R1.nDeployments_total, R1.det_ntrials(a_idx, 1, agent_db_idx, sensor_idx) / R1.nDeployments_total, upper(R1.detector_spread));

% MSE curves: {legend name, field, colour, line style}
C = orderedcolors("gem12");
methods = {
    'Proposed',    'avg_dep_mse',             C(1,:), '-o'
    'Undefended',  'avg_dep_mse_undefended',  C(4,:), '-s'
    'Oracle',      'avg_dep_mse_oracle',      C(5,:), '-^'
    'Full honest', 'avg_dep_mse_full_honest', C(7,:), ':'
    };
conv_names = string(fieldnames(R1.avg_dep_mse_conv)).';
conv_cols  = C([2 6 11 12], :);
for q = 1:numel(conv_names)
    methods(end+1, :) = {char(replace(conv_names(q), "_", " ")), "conv:" + conv_names(q), conv_cols(q,:), '--d'}; %#ok<SAGROW>
end
unknown = setdiff(show_methods, string(methods(:,1)));
if ~isempty(unknown)
    error('show_methods: not available in these results: %s (available: %s).', strjoin(unknown, ', '), strjoin(string(methods(:,1)), ', '));
end
methods = methods(ismember(string(methods(:,1)), show_methods), :);
mse_of = @(R, fld, prm, c) get_mse(R, fld, agent_db_idx, c, prm, scheme_idx, a_idx, sensor_idx) / norm_coeffs(prm);
across = @(fld, prm, c) cellfun(@(R) mse_of(R, fld, prm, c), L);   % 1 x n_lev

% Detection rates at the simulation's thresholds: {legend name, Pd(R, c), Pf(R, c), colour}
[~, ic]  = min(abs(R1.roc_corr_grid - R1.corr_min));
[~, ika] = min(abs(R1.roc_k_grid - R1.alpha_delta_k));
[~, ikt] = min(abs(R1.roc_k_grid - R1.t0_delta_k));
sp  = find(R1.roc_spread_names == R1.detector_spread, 1);
spd = replace(upper(R1.detector_spread), ["SN", "QN"], ["Sn", "Qn"]);
n_atk = @(R, c) R.roc_n_atk(a_idx, c, agent_db_idx, sensor_idx);
n_hon = @(R, c) R.roc_n_hon(a_idx, c, agent_db_idx, sensor_idx);
tests = {
    'Combined detector', @(R,c) R.det_tp(a_idx,c,agent_db_idx,sensor_idx) / (R.det_tp(a_idx,c,agent_db_idx,sensor_idx) + R.det_fn(a_idx,c,agent_db_idx,sensor_idx)), ...
                         @(R,c) R.det_fp(a_idx,c,agent_db_idx,sensor_idx) / n_hon(R,c), [0 0 0]
    sprintf('normcorr (corr_{min} = %g)', R1.corr_min), ...
                         @(R,c) R.roc_corr_det(ic,a_idx,c,agent_db_idx,sensor_idx) / n_atk(R,c), ...
                         @(R,c) R.roc_corr_fa(ic,a_idx,c,agent_db_idx,sensor_idx) / n_hon(R,c), C(2,:)
    sprintf('t_0 test (%s, k = %g)', spd, R1.t0_delta_k), ...
                         @(R,c) R.roc_t0_det(ikt,sp,a_idx,c,agent_db_idx,sensor_idx) / n_atk(R,c), ...
                         @(R,c) R.roc_t0_fa(ikt,sp,a_idx,c,agent_db_idx,sensor_idx) / n_hon(R,c), C(1,:)
    sprintf('\\alpha test (%s, k = %g)', spd, R1.alpha_delta_k), ...
                         @(R,c) R.roc_alpha_det(ika,sp,a_idx,c,agent_db_idx,sensor_idx) / n_atk(R,c), ...
                         @(R,c) R.roc_alpha_fa(ika,sp,a_idx,c,agent_db_idx,sensor_idx) / n_hon(R,c), C(3,:)
    };
rate_of = @(f, c) cellfun(@(R) f(R, c), L);   % 1 x n_lev
trial_rate = @(fld, R, c) R.(fld)(a_idx,c,agent_db_idx,sensor_idx) / R.det_ntrials(a_idx,c,agent_db_idx,sensor_idx);

% Near-oracle limit per channel SNR: worst (largest) CSI NMSE with proposed t0 within oracle_tol of oracle
near_oracle = nan(1, n_snr);
for c = snr_idx
    ok = across('avg_dep_mse', 2, c) <= (1 + oracle_tol) * across('avg_dep_mse_oracle', 2, c);
    last_ok = find(~ok, 1) - 1;          % levels are sorted from best to worst CSI
    if isempty(last_ok)
        last_ok = n_lev;
    end
    if last_ok >= 1
        near_oracle(c) = x_pos(last_ok);
    end
end

%% Figure 1: normalized MSE vs CSI NMSE (columns: t0, alpha; rows: channel SNR)
[tl, fig1] = grid_figure(1, numel(snr_idx), 2, {'Imperfect CSI: normalized MSE vs CSI NMSE', run_info});
cols = [2 1];   % t0 first
for c = snr_idx
    for r = 1:2
        prm = cols(r);
        ax = nexttile(tl);
        hold(ax, 'on');
        for q = 1:size(methods, 1)
            plot(ax, x_pos, across(methods{q,2}, prm, c), methods{q,4}, 'Color', methods{q,3}, ...
                'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', methods{q,1});
        end
        mark_limit(ax, near_oracle(c));
        ax.YScale = 'log';
        style_csi_axes(ax, x_pos, x_labels, font_size);
        if c == snr_idx(1)
            title(ax, sprintf('NMSE(%s)', param_lbl(prm)), 'FontWeight', 'normal');
        end
        if r == 1
            ylabel(ax, {sprintf('%g dB channel SNR', snr_vals(c)), 'normalized MSE'});
        end
        if c == snr_idx(end)
            xlabel(ax, 'CSI NMSE (dB)');
        end
        if r == leg_tile && c == snr_idx(1)
            legend(ax, 'Location', 'best', 'FontSize', font_size - 2);
        end
    end
end

%% Figure 2: Pd and Pf vs CSI NMSE (columns: Pd, Pf; rows: channel SNR)
[tl, fig2] = grid_figure(2, numel(snr_idx), 2, {'Imperfect CSI: detection vs CSI NMSE', run_info});
rate_lbl = ["P_d (per sensor)", "P_f (per sensor)"];
for c = snr_idx
    for col = 1:2   % 1: Pd, 2: Pf
        ax = nexttile(tl);
        hold(ax, 'on');
        for q = 1:size(tests, 1)
            plot(ax, x_pos, rate_of(tests{q, 1 + col}, c), '-o', 'Color', tests{q,4}, ...
                'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', tests{q,1});
        end
        mark_limit(ax, near_oracle(c));
        style_csi_axes(ax, x_pos, x_labels, font_size);
        if col == 2 && log_pf
            ax.YScale = 'log';
        else
            ylim(ax, [-0.05 1.08]);
            set_ticks(ax, 'y', 0:0.25:1);
        end
        if c == snr_idx(1)
            title(ax, rate_lbl(col), 'FontWeight', 'normal');
        end
        if col == 1
            ylabel(ax, {sprintf('%g dB channel SNR', snr_vals(c)), 'rate'});
        end
        if c == snr_idx(end)
            xlabel(ax, 'CSI NMSE (dB)');
        end
        if col == 2 && c == snr_idx(1)
            legend(ax, 'Location', 'best', 'FontSize', font_size - 2);
        end
    end
end

if save_figs
    experiment_name = "csi_sensitivity";
    file_info = attack + "_" + string(R1.agent_db_values(agent_db_idx)) + "_db_Tp=" + string(R1.Tp) + "_S=" + string(R1.sensor_vals(sensor_idx));
    save_figure(fig1, experiment_name, "mse_" + file_info);
    save_figure(fig2, experiment_name, "detection_" + file_info);
end

%% Table: near-oracle limit and the cost of CSI errors per channel SNR
fprintf('\nImperfect CSI (%s attack) | near-oracle: proposed t0 NMSE <= %.0f%% above the oracle\n', attack, 100*oracle_tol);
for c = snr_idx
    if isnan(near_oracle(c))
        lim = 'not reached';
    else
        lim = x_labels(x_pos == near_oracle(c)) + " dB";
    end
    fprintf('\n-- channel SNR %g dB | near-oracle down to CSI NMSE: %s --\n', snr_vals(c), lim);
    fprintf('%14s | %10s %10s %10s %10s | %9s %9s %9s %10s\n', 'CSI NMSE (dB)', 't0 prop', 't0/oracle', 'a prop', 'a/oracle', 'Pd', 'Pf', 'P(any FA)', 'P(fallbk)');
    t0p = across('avg_dep_mse', 2, c); t0o = across('avg_dep_mse_oracle', 2, c);
    ap  = across('avg_dep_mse', 1, c); ao  = across('avg_dep_mse_oracle', 1, c);
    pd  = rate_of(tests{1,2}, c);      pf  = rate_of(tests{1,3}, c);
    pfa = rate_of(@(R,cc) trial_rate('det_any_fa', R, cc), c);
    pfb = rate_of(@(R,cc) trial_rate('det_all_flagged', R, cc), c);
    for q = 1:n_lev
        fprintf('%14s | %10.3g %10.3g %10.3g %10.3g | %9.4f %9.4f %9.4f %10.5f\n', x_labels(q), t0p(q), t0p(q)/t0o(q), ap(q), ap(q)/ao(q), pd(q), pf(q), pfa(q), pfb(q));
    end
end


%% ---- Local functions ----
function v = get_mse(R, fld, ia, c, prm, is, ip, isen)
% Averaged MSE of one estimator; fld = "conv:<method>" selects a conventional baseline.
    if startsWith(fld, "conv:")
        A = R.avg_dep_mse_conv.(char(extractAfter(fld, "conv:")));
    else
        A = R.(char(fld));
    end
    v = A(1, ia, c, prm, is, ip, isen);
end

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

function style_csi_axes(ax, x_pos, x_labels, font_size)
    ax.FontSize = font_size;
    xlim(ax, [min(x_pos) - 2, max(x_pos) + 2]);
    ax.XTick = x_pos;
    ax.XTickLabel = x_labels;
    grid(ax, 'on');
    box(ax, 'on');
end

function mark_limit(ax, x_lim)
    if ~isnan(x_lim)
        xline(ax, x_lim, '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 1, 'HandleVisibility', 'off');
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
