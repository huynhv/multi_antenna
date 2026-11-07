% PLOT_ATTACKER_COUNT  Detection vs the number of compromised sensors (sensors 1..n).
% Run standalone, or after multi_antenna.m with experiment_list = "attacker_count" (detection only,
% one agent SNR and one channel SNR).
% Two figures per group of attacks in attack_groups (columns = the group's attacks):
%   Detection: rows Pd / Pf, x = number of attackers n. Pd and Pf per test at the simulation's
%             thresholds; dashed Pf = combined, no attack.
%   ROC:      rows = tests, one curve per n, 'o' at the current threshold.
% Attacks in hide_attacks are left out of every figure.
% A dotted vertical line marks n = S/2, where the robust median/spread estimates break down.

if exist('experiment', 'var') && experiment == "attacker_count"
    results_file = experiment + "_results.mat";
else
    results_file = "attacker_count_results.mat";
end
R = load(results_file);   % struct, so the caller's workspace isn't overwritten
experiment_name = erase(string(results_file), "_results.mat");   % used in the saved file names

font_size  = 12;     % font size of all axes (legends: 2 pt smaller)
save_figs  = true;   % save each figure as a 300 dpi PNG to ~\Desktop\Byzantine AirComp Results\<experiment>\
hide_attacks = ["peak"];   % attacks left out of every figure (still in the table)
attack_groups = {["noise", "sinusoid"], ["flip", "time shift"]};   % one detection + one ROC figure per group, columns = attacks

%% Sweep layout (single channel SNR / agent SNR / sensor count)
n_vals   = R.ac_counts(:).';
attack_groups = cellfun(@(g) g(~ismember(g, hide_attacks)), attack_groups, 'UniformOutput', false);
missing = setdiff([attack_groups{:}], R.ac_attacks);
if ~isempty(missing)
    error('attack_groups attack(s) not in these results: %s (attacks: %s).', strjoin(missing, ', '), strjoin(R.ac_attacks, ', '));
end
attack_groups = attack_groups(~cellfun(@isempty, attack_groups));
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

%% Detection figures (one per group in attack_groups): Pd and Pf vs number of attackers
row_lbl = ["P_d (per sensor)", "P_f (per sensor)"];
det_figs = cell(1, numel(attack_groups));
for g = 1:numel(attack_groups)
    atks = attack_groups{g};
    [tl, det_figs{g}] = grid_figure(g, 2, numel(atks), {'Detection vs number of compromised sensors', run_info});
    for row = 1:2
        for a = 1:numel(atks)
            ax = nexttile(tl);
            ax.FontSize = font_size;
            hold(ax, 'on');
            P = nan(numel(n_vals), 4);
            for q = 1:numel(n_vals)
                p = pivot_of(atks(a), n_vals(q));
                if row == 1
                    P(q,:) = pd_of(p);
                else
                    P(q,:) = pf_of(p);
                end
            end
            for k = 1:4
                plot(ax, n_vals, P(:,k), '-o', 'Color', test_colors(k,:), 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', test_names{k});
            end
            if row == 2
                pf_none = pf_of(p_none);
                yline(ax, pf_none(1), '--', 'Color', [0 0 0], 'DisplayName', 'Combined P_f, no attack');
            end
            xline(ax, S/2, ':', 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off');
            xlim(ax, [min(n_vals) - 0.3, max(n_vals) + 0.3]);
            ax.XTick = n_vals;
            ylim(ax, [-0.05 1.08]);
            set_ticks(ax, 'y', 0:0.25:1);
            grid(ax, 'on');
            box(ax, 'on');
            if row == 1
                title(ax, atks(a), 'FontWeight', 'normal');
            end
            if a == 1
                ylabel(ax, row_lbl(row));
            end
            if row == 2
                xlabel(ax, 'compromised sensors n');
            end
            if a == numel(atks)   % legend inside the rightmost plot of each row (each row has its own entries)
                legend(ax, 'Location', 'best', 'FontSize', font_size - 2);
            end
        end
    end
end

%% ROC figures (one per group in attack_groups): one curve per number of attackers
roc_tests = {
    'normcorr',                         @(p) R.roc_corr_det(:, p, 1, 1, 1) / n_atk(p),          @(p) R.roc_corr_fa(:, p, 1, 1, 1) / n_hon(p),          ic
    sprintf('\\alpha test (%s)', spd),  @(p) R.roc_alpha_det(:, sp, p, 1, 1, 1) / n_atk(p),     @(p) R.roc_alpha_fa(:, sp, p, 1, 1, 1) / n_hon(p),     ika
    sprintf('t_0 test (%s)', spd),      @(p) R.roc_t0_det(:, sp, p, 1, 1, 1) / n_atk(p),        @(p) R.roc_t0_fa(:, sp, p, 1, 1, 1) / n_hon(p),        ikt
    };
n_colors = parula(numel(n_vals) + 1);
n_tests  = size(roc_tests, 1);
roc_figs = cell(1, numel(attack_groups));
for g = 1:numel(attack_groups)
    atks = attack_groups{g};
    [tl, roc_figs{g}] = grid_figure(numel(attack_groups) + g, n_tests, numel(atks), ...
        {'ROC per test, one curve per number of compromised sensors (o: current threshold)', run_info});
    for k = 1:n_tests
        for a = 1:numel(atks)
            ax = nexttile(tl);
            ax.FontSize = font_size;
            hold(ax, 'on');
            for q = 1:numel(n_vals)
                p = pivot_of(atks(a), n_vals(q));
                pd = roc_tests{k,2}(p);
                pf = roc_tests{k,3}(p);
                plot(ax, pf, pd, '-', 'Color', n_colors(q,:), 'LineWidth', 1.5, 'DisplayName', sprintf('n = %d', n_vals(q)));
                i_now = roc_tests{k,4};
                plot(ax, pf(i_now), pd(i_now), 'o', 'Color', n_colors(q,:), 'MarkerFaceColor', n_colors(q,:), 'HandleVisibility', 'off');
            end
            plot(ax, [0 1], [0 1], ':', 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off');
            xlim(ax, [-0.05 1.05]);
            set_ticks(ax, 'x', 0:0.25:1);
            ylim(ax, [-0.05 1.08]);
            set_ticks(ax, 'y', 0:0.25:1);
            grid(ax, 'on');
            box(ax, 'on');
            if k == 1
                title(ax, atks(a), 'FontWeight', 'normal');
            end
            if a == 1
                ylabel(ax, {roc_tests{k,1}, 'P_d'});
            end
            if k == n_tests
                xlabel(ax, 'P_f');
            end
            if k == 1 && a == numel(atks)
                legend(ax, 'Location', 'southeast', 'FontSize', font_size - 2);   % ROC curves leave the lower right empty
            end
        end
    end
end

if save_figs
    file_info = string(R.ac_agent_db) + "_db_channel_" + string(R.ac_channel_snr) + "_db_Tp=" + string(R.Tp) + "_S=" + string(S);
    for g = 1:numel(attack_groups)
        group_tag = strjoin(replace(attack_groups{g}, " ", "_"), "_");
        save_figure(det_figs{g}, experiment_name, "detection_" + group_tag + "_" + file_info);
        save_figure(roc_figs{g}, experiment_name, "roc_" + group_tag + "_" + file_info);
    end
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
