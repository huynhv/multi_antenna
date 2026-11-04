% PLOT_ATTACKER_COUNT  Detection vs the number of compromised sensors (sensors 1..n).
% Run standalone, or after multi_antenna.m with experiment_list = "attacker_count" (detection only,
% one agent SNR and one channel SNR).
%   Figure 1: rows Pd / Pf / per-trial rates, columns attack types, x = number of attackers n.
%             Pd and Pf per test at the simulation's thresholds; dashed Pf = combined, no attack;
%             per-trial: P(every attacker flagged) and P(at least one honest sensor flagged).
%   Figure 2: ROC per test (columns) and attack (rows), one curve per n, 'o' at the current threshold.
% A dotted vertical line marks n = S/2, where the robust median/spread estimates break down.

if exist('experiment', 'var') && experiment == "attacker_count"
    results_file = experiment + "_results.mat";
else
    results_file = "attacker_count_results.mat";
end
R = load(results_file);   % struct, so the caller's workspace isn't overwritten
experiment_name = erase(string(results_file), "_results.mat");   % used in the saved file names

save_figs = true;   % save each figure as a 300 dpi PNG to ~\Desktop\Byzantine AirComp Results\<experiment>\

%% Sweep layout (single channel SNR / agent SNR / sensor count)
n_vals   = R.ac_counts(:).';
att_list = R.ac_attacks;
p_none   = find(R.pivot_types == "none", 1);
pivot_of = @(a, n) find(R.pivot_types == a & R.pivot_n == n, 1);   % sweep entry of (attack, n)
S        = R.sensor_vals(1);
C        = orderedcolors("gem12");
run_info = sprintf('sensors 1..n compromised of S = %d | agent SNR %g dB, channel SNR %g dB | detector: normcorr + %s (k = %g / %g), corr_{min} = %g', ...
    S, R.ac_agent_db, R.ac_channel_snr, upper(R.detector_spread), R.alpha_delta_k, R.t0_delta_k, R.corr_min);

[~, ic]  = min(abs(R.roc_corr_grid - R.corr_min));
[~, ika] = min(abs(R.roc_k_grid - R.alpha_delta_k));
[~, ikt] = min(abs(R.roc_k_grid - R.t0_delta_k));
sp = find(R.roc_spread_names == R.detector_spread, 1);
spd = replace(upper(R.detector_spread), ["SN", "QN"], ["Sn", "Qn"]);

% Rates of one sweep entry p: Pd and Pf for {combined, normcorr, t0 test, alpha test}
n_atk = @(p) R.roc_n_atk(p, 1, 1, 1);
n_hon = @(p) R.roc_n_hon(p, 1, 1, 1);
pd_of = @(p) [R.det_tp(p,1,1,1) / (R.det_tp(p,1,1,1) + R.det_fn(p,1,1,1)), ...
              R.roc_corr_det(ic, p, 1, 1, 1) / n_atk(p), ...
              R.roc_t0_det(ikt, sp, p, 1, 1, 1) / n_atk(p), ...
              R.roc_alpha_det(ika, sp, p, 1, 1, 1) / n_atk(p)];
pf_of = @(p) [R.det_fp(p,1,1,1) / n_hon(p), ...
              R.roc_corr_fa(ic, p, 1, 1, 1) / n_hon(p), ...
              R.roc_t0_fa(ikt, sp, p, 1, 1, 1) / n_hon(p), ...
              R.roc_alpha_fa(ika, sp, p, 1, 1, 1) / n_hon(p)];
test_names  = {'Combined detector', 'normcorr', sprintf('t_0 test (%s)', spd), sprintf('\\alpha test (%s)', spd)};
test_colors = [0 0 0; C(2,:); C(1,:); C(3,:)];

%% Figure 1: Pd, Pf and per-trial rates vs number of attackers
[tl, fig1] = grid_figure(1, 3, numel(att_list), {'Detection vs number of compromised sensors', run_info});
row_lbl = ["P_d (per sensor)", "P_f (per sensor)", "per-trial rate"];
for row = 1:3
    for a = 1:numel(att_list)
        ax = nexttile(tl);
        hold(ax, 'on');
        P = nan(numel(n_vals), 4);
        trial_rates = nan(numel(n_vals), 2);
        for q = 1:numel(n_vals)
            p = pivot_of(att_list(a), n_vals(q));
            if row == 1
                P(q,:) = pd_of(p);
            elseif row == 2
                P(q,:) = pf_of(p);
            else
                trial_rates(q,:) = [R.det_all_caught(p,1,1,1), R.det_any_fa(p,1,1,1)] / R.det_ntrials(p,1,1,1);
            end
        end
        if row < 3
            for k = 1:4
                plot(ax, n_vals, P(:,k), '-o', 'Color', test_colors(k,:), 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', test_names{k});
            end
            if row == 2
                pf_none = pf_of(p_none);
                yline(ax, pf_none(1), '--', 'Color', [0 0 0], 'DisplayName', 'Combined P_f, no attack');
            end
        else
            plot(ax, n_vals, trial_rates(:,1), '-o', 'Color', [0 0 0], 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', 'P(all attackers flagged)');
            plot(ax, n_vals, trial_rates(:,2), '--s', 'Color', C(7,:), 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', 'P(\geq 1 honest flagged)');
        end
        xline(ax, S/2, ':', 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off');
        xlim(ax, [min(n_vals) - 0.3, max(n_vals) + 0.3]);
        ax.XTick = n_vals;
        ylim(ax, [-0.05 1.08]);
        grid(ax, 'on');
        box(ax, 'on');
        if row == 1
            title(ax, att_list(a), 'FontWeight', 'normal');
        end
        if a == 1
            ylabel(ax, row_lbl(row));
        end
        if row == 3
            xlabel(ax, 'compromised sensors n');
        end
        if a == numel(att_list)
            legend(ax, 'Location', 'eastoutside', 'FontSize', 9);
        end
    end
end

%% Figure 2: ROC per test, one curve per number of attackers
roc_tests = {
    'normcorr',                         @(p) R.roc_corr_det(:, p, 1, 1, 1) / n_atk(p),          @(p) R.roc_corr_fa(:, p, 1, 1, 1) / n_hon(p),          ic
    sprintf('\\alpha test (%s)', spd),  @(p) R.roc_alpha_det(:, sp, p, 1, 1, 1) / n_atk(p),     @(p) R.roc_alpha_fa(:, sp, p, 1, 1, 1) / n_hon(p),     ika
    sprintf('t_0 test (%s)', spd),      @(p) R.roc_t0_det(:, sp, p, 1, 1, 1) / n_atk(p),        @(p) R.roc_t0_fa(:, sp, p, 1, 1, 1) / n_hon(p),        ikt
    };
n_colors = parula(numel(n_vals) + 1);
[tl, fig2] = grid_figure(2, numel(att_list), size(roc_tests, 1), {'ROC per test, one curve per number of compromised sensors (o: current threshold)', run_info});
for a = 1:numel(att_list)
    for k = 1:size(roc_tests, 1)
        ax = nexttile(tl);
        hold(ax, 'on');
        for q = 1:numel(n_vals)
            p = pivot_of(att_list(a), n_vals(q));
            pd = roc_tests{k,2}(p);
            pf = roc_tests{k,3}(p);
            plot(ax, pf, pd, '-', 'Color', n_colors(q,:), 'LineWidth', 1.5, 'DisplayName', sprintf('n = %d', n_vals(q)));
            i_now = roc_tests{k,4};
            plot(ax, pf(i_now), pd(i_now), 'o', 'Color', n_colors(q,:), 'MarkerFaceColor', n_colors(q,:), 'HandleVisibility', 'off');
        end
        plot(ax, [0 1], [0 1], ':', 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off');
        xlim(ax, [-0.05 1.05]);
        ylim(ax, [-0.05 1.08]);
        grid(ax, 'on');
        box(ax, 'on');
        if a == 1
            title(ax, roc_tests{k,1}, 'FontWeight', 'normal');
        end
        if k == 1
            ylabel(ax, {char(att_list(a)), 'P_d'});
        end
        if a == numel(att_list)
            xlabel(ax, 'P_f');
        end
        if a == 1 && k == size(roc_tests, 1)
            legend(ax, 'Location', 'eastoutside', 'FontSize', 9);
        end
    end
end

if save_figs
    file_info = string(R.ac_agent_db) + "_db_channel_" + string(R.ac_channel_snr) + "_db_Tp=" + string(R.Tp) + "_S=" + string(S);
    save_figure(fig1, experiment_name, "detection_" + file_info);
    save_figure(fig2, experiment_name, "roc_" + file_info);
end

%% Table
fprintf('\nAttacker count (detection only): agent SNR %g dB, channel SNR %g dB, S = %d\n', R.ac_agent_db, R.ac_channel_snr, S);
disp(R.summary_tbl)


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
    w = min(scr(3) - 80,  260*cols + 260);
    h = min(scr(4) - 140, 230*rows + 140);
    set(fig, 'Position', [40 + 25*fig_count, max(40, scr(4) - h - 90 - 25*fig_count), w, h]);
    tl = tiledlayout(fig, rows, cols, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, ttl{1}, 'FontWeight', 'bold');
    subtitle(tl, ttl{2});
end
