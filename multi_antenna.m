%%% 10-6-2026 --> revisit tomorrow

% Clear all variables and close all existing figures.
clearvars
close all

addpath('.\functions')
verbosity_level = 5;   % 1 = major sections only, 2 = + sub-steps/channel-SNR, 3 = + diagnostics

% Set seed for reproducibility
seed = 264;
rng(seed);

use_W = true;

nDeployments_total = 100;   % total deployments to simulate
dep_batch_size     = 10; % 10;    % deployments held in memory at once (~3-4 GB per batch of 10)
% Set the number of measurements per sensor, starting from 0
K = 301;
% Set number of trials.
nTrials = 500;

% Set number of sensors
sensor_vals = 9;
sensor_dimension = length(sensor_vals);

% Set observation window length in seconds.
T0 = 4;

% Define event signal parameters.
Tp = 1;

% Define time domain parameters
dt = T0/(K-1);
t = dt*(0:K-1)';

offset_idx = ceil(Tp/dt);

% Define frequency domain parameters
fs = 1/dt;
L = K-1;
f = 0:fs/(2*L):fs/2;
omega = 2*pi*f;

% Calculate the signal energy (such that the signal is assumed to be within the
% interior to the observation period).
func = @(t) abs(sensor_signal(t, Tp, 1)).^2;
Ps = integral(func,0,T0)/T0;

% Define true parameter values.
alpha_true = 3;
t0_true = dt*ceil(1.5/dt);

% Estimation mode: false = joint (alpha evaluated at the estimated t0),
%                  true  = alpha evaluated at the true t0.
alpha_uses_true_t0 = false;

% Set the number of antennas
num_antennas = 9;
Mtot = 2*num_antennas;

%%% Always generate w, n, and channel gain for maximum S, then use subsets for different S trials
S_max = max(sensor_vals);
% Channel gains for every deployment (small). The large noise arrays are drawn per batch.
all_gi_total = (randn(1,S_max,num_antennas,1,nDeployments_total) + 1i*randn(1,S_max,num_antennas,1,nDeployments_total))/sqrt(2);

%% Generate mi and ti
max_ti = 1;
min_ti = 0.8; 
path_loss_exp = 0.5;
shadow_sigma_dB = 0.1;

m_ref = 1;              % reference amplitude at distance min_ti
c_gain = m_ref * min_ti^path_loss_exp;   % anchors: unshadowed mi = m_ref at d = min_ti

all_di = sort(unifrnd(min_ti, max_ti, 1, S_max, 1, 1, nDeployments_total), 2, 'ascend');

shadowing = 10.^(shadow_sigma_dB/20 * randn(1,S_max,1,1,nDeployments_total));
all_mi_total = c_gain ./ (all_di.^path_loss_exp) .* shadowing;

%%
selected_schemes = "EPC";

% Define Rayleigh distribution parameters.
rayleigh_factor = 1/sqrt(2);
E_mag_g_sqr = 2*rayleigh_factor^2;
E_mi_sqr = mean(all_mi_total(:).^2); % --> could generate empirical value here

% Set miscellaneous parameters
norm_fact = 1;
epsilon = 1e-4;
Bee = pi/Tp;

% Set maximum value of t0 such that the signal from each sensor is guaranteed to
% be in the interior of the observation period.
t0_max = T0-Tp-max_ti;

if t0_true > t0_max
    error("True value of t0 exceeds maximum!")
end

% Keep event timing on the sample grid (avoids a t0-MSE floor of ~(dt/2)^2
% and interpolated time-shift attacks).
on_grid = @(x) abs(x/dt - round(x/dt)) < 1e-9;
assert(on_grid(t0_true), 't0_true = %g is not a multiple of dt = %g', t0_true, dt);
assert(on_grid(Tp),      'Tp = %g is not a multiple of dt = %g', Tp, dt);

%% Byzantine attacker configuration
attacker_enabled = true;
attacker_idx     = [1,2,3]; % [2,4]      % sensor index/indices (within 1:S) that are compromised
attacker_db      = 0;      % attacker "SNR", same convention as agent_db (see below)
attacker_seed = 42;
% Attack-specific parameters (each is used only by its attack type)
attacker_beta = 0.1;   % "amplitude scale": attacker transmits beta * u_i(t)
t0_attacker   = 3;     % "time shift": attacker's pulse is moved to start at t0_attacker

% Time-shift sensitivity experiment (experiment_list = "time_shift_sensitivity"): one time-shift attack
% per multiple n, with the attacker's pulse moved to t0_true + n*ts_delta. A no-attack reference
% ("none") is added automatically. Fractional n gives sub-sample shifts (linear interpolation of the
% smooth matched-filter output in shift_waveform).
ts_delta     = dt;                                % Delta = one sample
ts_multiples = [0.25 0.5 1 2 4 8 16 32 64];       % n

% Spread estimator used by the per-sensor alpha/t0 outlier tests IN THE DETECTOR (all three are
% always computed for the ROC): "mad" (1.4826*MAD), "sn" or "qn" (Rousseeuw-Croux).
detector_spread = "qn";
assert(ismember(detector_spread, ["mad", "sn", "qn"]), 'detector_spread must be "mad", "sn" or "qn".');

% Conventional baselines built from the per-sensor estimates (alpha_i_hat, t0_i_hat), independent of
% the proposed detector. Any subset of:
%   Simple robust aggregation (coordinate-wise, no outlier tests):
%     "median"        : median over all sensors
%     "trimmed_mean"  : mean after dropping the conv_trim smallest and conv_trim largest values
%   Optional ablations (MAD outlier rejection, then averaging):
%     "skipped_mean"  : mean of the sensors kept by a conventional MAD outlier rule (median if none kept)
%     "skipped_wmean" : inverse-variance weighted mean of the same kept sensors (same fallback)
conventional_methods = ["median"];
valid_conventional = ["median", "trimmed_mean", "skipped_mean", "skipped_wmean"];
assert(all(ismember(conventional_methods, valid_conventional)), ...
    'conventional_methods may only contain: %s.', strjoin(valid_conventional, ', '));
conventional_methods = unique(conventional_methods, 'stable');
% Trimmed mean: values trimmed from EACH end. Default f = number of compromised sensors (the textbook
% Byzantine choice and the trimmed mean's best case -- it assumes the attack size is known).
conv_trim = numel(attacker_idx);
% Rejection rule of the skipped means: drop sensor i if |x_i - median| / (1.4826*MAD) > conv_mad_k
%   "per_parameter": alpha outliers dropped from the alpha mean, t0 outliers from the t0 mean (Hampel X84)
%   "joint"        : drop the sensor from both means if either its alpha or its t0 is an outlier
conv_mad_k     = 3;
conv_rejection = "per_parameter";
assert(ismember(conv_rejection, ["per_parameter", "joint"]), 'conv_rejection must be "per_parameter" or "joint".');

w_psd_constant = 1;
n_psd_constant = 1;
N = 2*K - 1;

%% Subspace-projection residual test: honest-signal basis Phi and its
t0_grid_top = reshape(t,1,1,[]);
[~, R00_full] = mf_integral_fft(sensor_signal(t-t0_grid_top,Tp,norm_fact), sensor_signal(t,Tp,norm_fact), 1, 1, K, dt, Tp);
D_template = reshape(R00_full(:,:,:,1), K, K);   % K x K, [D_template]_{k,l} = rho(t_k - t_l)

[U_D, Sigma_D, ~] = svd(D_template, 'econ');
% sv_energy = cumsum(diag(Sigma_D).^2) / sum(diag(Sigma_D).^2);
% d_sub = find(sv_energy >= 0.999, 1, 'first');
% Phi = U_D(:, 1:d_sub);           % K x d_sub: honest-signal subspace

% Eigenvalues of D_template itself (already computed via U_D/Sigma_D above)
lam_D = diag(Sigma_D);   % K x 1
% Project D_template's columns into U_D's eigenbasis once -- reused every
% (sensor, deployment) below since it doesn't depend on either.
D_template_v = U_D.' * D_template;   % K x K
% d rho(t - tau)/d tau for every candidate tau (central differences across D_template's columns),
% projected the same way -- gives the per-sensor t0_i precision used by the weighted skipped mean.
D_deriv = zeros(K, K);
D_deriv(:, 2:K-1) = (D_template(:, 3:K) - D_template(:, 1:K-2)) / (2*dt);
D_deriv(:, 1)     = (D_template(:, 2) - D_template(:, 1)) / dt;
D_deriv(:, K)     = (D_template(:, K) - D_template(:, K-1)) / dt;
D_deriv_v = U_D.' * D_deriv;          % K x K
mfTemplateFFT_raw = fft(sensor_signal(t, Tp, norm_fact), N, 1);   % used by estimate_with_geometry

% log_msg(verbosity_level, 1, 'Subspace-projection residual test: d = %d, K-d = %d', d_sub, K-d_sub);

%% Define Anonymous Expressions
mag_sqr_S0_internal = @(w) (norm_fact^2) * (2*Bee.^2 .* (1 + cos(w.*pi./Bee)) ) ./ (w.^2 - Bee.^2).^2;
mag_sqr_S0 = @(w) zeroIfnan(mag_sqr_S0_internal(w)) + (w == Bee | w == -Bee) .* (mag_sqr_S0_internal(w-epsilon) + mag_sqr_S0_internal(w+epsilon))./2;

%%% Expressions for multi-antenna
Vm_arr_internal = @(w, lambda) lambda .* mag_sqr_S0_internal(w) + 1;
Vm_arr = @(w, lambda) zeroIfnan(Vm_arr_internal(w, lambda)) + (w == Bee | w == -Bee) .* (Vm_arr_internal(w-epsilon, lambda) + Vm_arr_internal(w+epsilon, lambda))./2;

Hm_arr_internal = @(w, lambda) (1 + epsilon) ./ (Vm_arr_internal(w, lambda) + epsilon);
Hm_arr = @(w, lambda) zeroIfnan(Hm_arr_internal(w, lambda)) + (w == Bee | w == -Bee) .* (Hm_arr_internal(w-epsilon, lambda) + Hm_arr_internal(w+epsilon, lambda))./2;

% alpha_branch_fisher = @(bm, lambda) integral(@(w) bm.^2 .* mag_sqr_S0(w).^2 ./ Vm_arr(w, lambda), -Inf, Inf, 'ArrayValued', true);
% t0_branch_fisher = @(bm, lambda) integral(@(w) bm.^2 .* (w .* mag_sqr_S0(w)).^2 ./ Vm_arr(w, lambda), -Inf, Inf, 'ArrayValued', true);

%% Define experiment list
experiment_list = "time_shift_sensitivity";   % and/or "multi_antenna_base", "time_shift_sensitivity"

for experiment_idx = 1:numel(experiment_list)
    experiment = experiment_list(experiment_idx);

    pad_str = '----------';
    msg = [pad_str  ' ' 'Running: ' char(experiment) ' ' pad_str];
    log_msg(verbosity_level, 1, '%s', msg);

    % Configure strategies
    agent_db_values = 0;
    channel_snr = linspace(0,10,5);

    % Attacks to test -- each runs the full Monte Carlo sweep. Valid names are listed below.
    attacks = ["none", "noise", "sinusoid", "peak", "flip", "time shift"]; % ["none", "noise", "sinusoid", "peak", "flip", "time shift"];

    valid_attacks = ["none", "noise", "sinusoid", "peak", "flip", "time shift", "amplitude scale"];
    attacks = unique(attacks, 'stable');   % no duplicate runs
    if isempty(attacks)
        error("No attack type has been specified!")
    end
    bad = attacks(~ismember(attacks, valid_attacks));
    if ~isempty(bad)
        error("Unknown attack type(s): %s. Valid: %s", strjoin(bad, ", "), strjoin(valid_attacks, ", "));
    end

    % Each sweep entry ("pivot") has a label (pivot_vals, used in tables and plots), an attack type
    % (pivot_types) and, for time shifts, the shift in seconds (pivot_shift).
    if experiment == "time_shift_sensitivity"
        ts_shift_s   = ts_multiples * ts_delta;          % shift of each attack (s)
        ts_shift_pct = 100 * ts_shift_s / t0_true;       % ... as a percentage of t0_true
        attacks      = ["none", "time shift " + compose("%g", ts_multiples) + "dt"];
        pivot_types  = ["none", repmat("time shift", 1, numel(ts_multiples))];
        pivot_shift  = [NaN, ts_shift_s];
        log_msg(verbosity_level, 1, 'Time-shift sensitivity: Delta = dt = %.4g ms = %.4g%% of t0_true (%.4g s)', ...
            1e3*ts_delta, 100*ts_delta/t0_true, t0_true);
        for q = 1:numel(ts_multiples)
            log_msg(verbosity_level, 1, '  n = %-5g shift = %8.3f ms = %7.3f%% of t0_true', ts_multiples(q), 1e3*ts_shift_s(q), ts_shift_pct(q));
        end
    else
        pivot_types = attacks;
        pivot_shift = nan(size(attacks));
        pivot_shift(attacks == "time shift") = t0_attacker - t0_true;
    end
    pivot_vals = attacks;
    
    constraints = "none";
    transforms = "none";
    coeff_types = "none"; 
    approaches = "none";
    
    iter_arr = [];
    all_strats = [];

    for i = 1:numel(constraints)
        for j = 1:numel(approaches)
            for k = 1:numel(coeff_types)
                for l = 1:numel(transforms)
                    if i == 1
                        iter_arr = [iter_arr; coeff_types(k) + " " + transforms(l) + " " + approaches(j)];
                    end
                    all_strats = [all_strats; constraints(i) + " " + coeff_types(k) + " " + transforms(l) + " " + approaches(j)];
                end
            end
        end
    end
    
    iter_length = length(transforms) * length(coeff_types) * length(approaches);
    
    num_strats = length(all_strats);
    
    % Duplicate channel snr matrices for experiments.
    channel_db_values = repmat(channel_snr, length(selected_schemes), 1, length(agent_db_values));
    
    % Initialize empty performance metric matrices
    rho_crlb = zeros(num_strats,numel(agent_db_values),size(channel_db_values,2),2,numel(selected_schemes),numel(pivot_vals),sensor_dimension,nDeployments_total);
    rho_crlb_oracle = rho_crlb;   % Theoretical CRLB with the attacker's channel treated as absent -- upper bound, only under attack
    
    rho_empirical_bias = rho_crlb;
    rho_empirical_var = rho_crlb;
    rho_empirical_mse = rho_crlb;
    rho_empirical_mse_full_honest = rho_crlb;   % Full S-sensor array, attacker never hijacked anyone -- "what if there had been no attack at all"
    rho_empirical_mse_undefended = rho_crlb;   % MSE with NO defense applied -- only meaningful/populated under attack
    rho_empirical_mse_oracle = rho_crlb;   % MSE with attacker PERFECTLY, freely removed -- upper bound, only under attack
    rho_empirical_mse_basis_supp = rho_crlb;
    rho_mse_conv = struct();   % one MSE array per conventional baseline (same layout as rho_empirical_mse)
    for cm = conventional_methods
        rho_mse_conv.(char(cm)) = rho_crlb;
    end

    % Detection statistics per (attack, channel SNR, agent SNR, sensor count), summed over trials and deployments.
    det_size = [numel(pivot_vals), size(channel_db_values,2), length(agent_db_values), sensor_dimension];
    det_tp = zeros(det_size);  det_fp = zeros(det_size);  det_fn = zeros(det_size);
    fp_rate_combined = zeros(det_size);   % honest-sensor false-flag rate

    % ROC (Pd vs Pf): per-sensor flag counts with each test's OWN threshold swept over its
    % natural range (no normalization). The +/-Inf end points flag every sensor / none, so
    % each curve can be closed at (1,1) and (0,0).
    %   normcorr: flagged if norm_corr < corr_min                              -> sweep corr_min
    %   alpha/t0: flagged if |x_i - median| / scale > k, scale = 1.4826*MAD, Sn or Qn  -> sweep k
    roc_corr_grid = [-Inf, linspace(0, 1, 201), Inf];   % corr_min values
    roc_k_grid    = [-Inf, 0:0.125:10, Inf];            % threshold k in robust-sigma units (alpha_delta_k / t0_delta_k)
    roc_spread_names = ["mad", "sn", "qn"];             % 2nd dimension of the alpha/t0 arrays
    roc_corr_det  = zeros([numel(roc_corr_grid), det_size]);   % attacker cases flagged
    roc_corr_fa   = zeros([numel(roc_corr_grid), det_size]);   % honest cases flagged
    roc_alpha_det = zeros([numel(roc_k_grid), numel(roc_spread_names), det_size]);
    roc_alpha_fa  = zeros([numel(roc_k_grid), numel(roc_spread_names), det_size]);
    roc_t0_det    = zeros([numel(roc_k_grid), numel(roc_spread_names), det_size]);
    roc_t0_fa     = zeros([numel(roc_k_grid), numel(roc_spread_names), det_size]);
    roc_n_atk = zeros(det_size);   % number of attacker-sensor cases (sensor x trial x deployment)
    roc_n_hon = zeros(det_size);   % number of honest-sensor cases

    % Start run timer
    log_msg(verbosity_level, 1, 'Starting runtime...');
    loopTic = tic;

    nBatches = ceil(nDeployments_total / dep_batch_size);
    batch_time = nan(1, nBatches);   % wall-clock seconds per deployment batch
    for batch_idx = 1:nBatches
        batch_start = tic;
        dep_range = (batch_idx-1)*dep_batch_size + 1 : min(batch_idx*dep_batch_size, nDeployments_total);
        nDeployments = numel(dep_range);   % deployments in THIS batch -- every array below is sized by it
        log_msg(verbosity_level, 1, 'Deployment batch %d of %d (deployments %d-%d)', batch_idx, nBatches, dep_range(1), dep_range(end));

        % Small per-deployment quantities are sliced; large noise arrays are drawn per batch.
        all_gi = all_gi_total(:,:,:,:,dep_range);
        all_mi = all_mi_total(:,:,:,:,dep_range);
        rng(seed + batch_idx);
        all_w = randn(K,S_max,1,nTrials,nDeployments);                     % sensor noise
        n = randn(K,1,num_antennas,nTrials,nDeployments) + 1i*randn(K,1,num_antennas,nTrials,nDeployments);   % server noise

        %% Build every attack's waveform once (depends on the attack only). Each attack
        % re-seeds with attacker_seed, so its draws don't depend on which other attacks
        % are listed or in what order.
        f_nominal = Bee/(2*pi);   % = 1/(2*Tp), reference frequency near the honest passband
        attack_waveforms = cell(numel(pivot_vals), 1);   % unit-energy K x S_max x 1 x nTrials x D, or [] if unused
        for pivot_idx = 1:numel(pivot_vals)
            if pivot_types(pivot_idx) == "noise"
                rng(attacker_seed + batch_idx);   % fresh attack noise for each batch of deployments
            else
                rng(attacker_seed);               % same deterministic waveform in every batch
            end
            switch pivot_types(pivot_idx)
                case "noise"
                    wf = randn(K,S_max,1,nTrials,nDeployments);

                case "sinusoid"
                    f_c   = f_nominal * (0.5 + rand(1,S_max));
                    phase = 2*pi*rand(1,S_max);
                    wf = repmat(sin(2*pi*t.*f_c + phase), 1,1,1,nTrials,nDeployments);   % K x S_max, repeated

                case "peak"
                    peak_same_location = false;
                    candidate_locs = t(2:end);   % the two worst-case candidates identified earlier
                    if peak_same_location
                        peak_loc = repmat(candidate_locs(randi(numel(candidate_locs))), 1, S_max);
                    else
                        peak_loc = candidate_locs(randi(numel(candidate_locs), 1, S_max));
                    end
                    waveform = zeros(K, S_max);
                    for a = 1:S_max
                        [~, idx] = min(abs(t - peak_loc(a)));
                        waveform(idx, a) = 1;
                    end
                    wf = repmat(waveform, 1,1,1,nTrials,nDeployments);

                otherwise   % "flip", "amplitude scale", "time shift" transform the honest ui instead
                    wf = [];
            end
            if ~isempty(wf)   % normalize to unit total (dt-weighted) energy per slot
                wf = wf ./ sqrt(sum(wf(:,:,1,1,1).^2, 1) * dt);
            end
            attack_waveforms{pivot_idx} = wf;
        end
        rng(seed);

        % Attacker amplitude, same convention as the sensor noise with attacker_db in place of agent_db.
        gamma_w_attacker = (dt/w_psd_constant) * (Ps / db2magTen(attacker_db));
        attacker_scale = sqrt(gamma_w_attacker * w_psd_constant / dt);

        ai = 1;
        bi = 0;

        for agent_db_idx = 1:length(agent_db_values)
            agent_db = agent_db_values(agent_db_idx);

            % Set sensor noise power. Note the factor of 1/dt, which is equivalent to
            % applying an anti-aliasing filter.
            gamma_w = (dt/w_psd_constant) * (Ps / db2magTen(agent_db));
            scaled_w = sqrt(gamma_w * w_psd_constant /dt) * all_w;   % variance gamma_w/dt
            scaled_w_psd_constant = gamma_w * w_psd_constant;       % = Nw/2

            for sensor_idx = 1:length(sensor_vals)
                S = sensor_vals(sensor_idx);
                save_dim = sensor_idx;

                % Select subset of mi and gi values.
                mi_5d = all_mi(:,1:S,:,:,:);
                mi_4d = reshape(mi_5d, 1, S, 1, nDeployments);
                gi = all_gi(:,1:S,:,:,:);

                % Honest matched-filter outputs (attack-independent -- computed once).
                noisy_unified_xi = alpha_true .* (mi_5d .* sensor_signal(t-t0_true, Tp, norm_fact)) + scaled_w(:,1:S,:,:,:);
                [~, ui_no_attack] = mf_integral_fft(noisy_unified_xi, mi_5d .* sensor_signal(t, Tp, norm_fact), 1, 1, K, dt, Tp);

                %% What each attack makes the compromised sensors transmit in place of their honest ui
                atk = attacker_idx(attacker_idx <= S);
                ui_honest_atk = ui_no_attack(:,atk,:,:,:);
                ui_atk_by_attack = cell(numel(pivot_vals), 1);
                for pivot_idx = 1:numel(pivot_vals)
                    switch pivot_types(pivot_idx)
                        case "flip"              % -alpha*m_i^2*rho(t-t0) - w_i(t)
                            ui_atk_by_attack{pivot_idx} = -ui_honest_atk;
                        case "amplitude scale"   % beta*(alpha*m_i^2*rho(t-t0) + w_i(t))
                            ui_atk_by_attack{pivot_idx} = attacker_beta * ui_honest_atk;
                        case "time shift"        % shifted copy, truncated to [0, T0]
                            ui_atk_by_attack{pivot_idx} = shift_waveform(ui_honest_atk, t, pivot_shift(pivot_idx));
                        case "none"              % no attack: compromised sensors stay honest
                            ui_atk_by_attack{pivot_idx} = [];   % unused -- out_cws is the honest signal
                        otherwise                % replaced by a scaled attack waveform
                            ui_atk_by_attack{pivot_idx} = attacker_scale * attack_waveforms{pivot_idx}(:,atk,:,:,:);
                    end
                end

                %% Define search space for t0 (single-antenna baseline only)
                if ~use_W
                    t0_search_space = repmat(reshape(t,1,1,[]), 1,1,1,nDeployments);
                    [~, ui_t0_search] = mf_integral_fft(mi_4d .* sensor_signal(t-t0_search_space, Tp, norm_fact), mi_4d .* sensor_signal(t, Tp, norm_fact), 1, 1, K, dt, Tp);
                    mfTemplateFFT = fft(mi_4d .* sensor_signal(t, Tp, norm_fact), N, 1);
                end

                for scheme_idx = 1:length(selected_schemes)
                    scheme = selected_schemes(scheme_idx);

                    if scheme == "EPC"
                        % Apply exact phase compensation
                        ref_ant_idx = 1;
                        g_ref = gi(:, :, ref_ant_idx, :, :);
                        phase_comp = exp(-1j * angle(g_ref));
                        g_tilde = gi .* phase_comp;
                    end

                    % Sensor combining at each antenna: one (K x s)*(s x M) product per (trial, deployment)
                    % page, instead of broadcasting a K x S x M x nTrials x D temporary (~2 GB).
                    combine = @(u, g_sub) reshape(pagemtimes(reshape(u, K, size(u,2), nTrials, nDeployments), ...
                        reshape(g_sub, size(u,2), num_antennas, 1, nDeployments)), K, 1, num_antennas, nTrials, nDeployments);
                    out_cws_no_attack  = combine(ui_no_attack, g_tilde);
                    out_cws_honest_atk = combine(ui_honest_atk, g_tilde(:,atk,:,:,:));   % compromised sensors' honest share

                    %% Channel-independent whitening basis (eigenvalues rescaled per channel SNR below)
                    m = reshape(mi_5d, S, nDeployments);                   % S x nDeployments
                    g = reshape(g_tilde, S, num_antennas, nDeployments);
                    G_R = permute(cat(2, real(g), imag(g)), [2,1,3]);      % Mtot x S x nDeployments
                    [UB_transpose, lambda_B_base] = whitening_basis(G_R, m, gamma_w);

                    % Oracle CRLB: same construction using only the honest sensors.
                    if attacker_enabled
                        honest_only = setdiff(1:S, attacker_idx);
                        S_honest   = numel(honest_only);
                        m_oracle   = m(honest_only,:);
                        G_R_oracle = G_R(:, honest_only, :);
                        [UB_transpose_oracle, lambda_B_base_oracle] = whitening_basis(G_R_oracle, m_oracle, gamma_w);
                        mu_oracle = m_oracle.^2;
                    end

                    %% --- Precompute null-space isolation projectors (channel-independent;
                    % reused across every channel-SNR point, attack and trial for this S/scheme) ---
                    r_dim = Mtot - S + 1;   % surviving real dimensions after nulling S-1 sensors
                    Pi_all = cell(S,1);
                    breve_g_all = cell(S,1);
                    for i = 1:S
                        others = setdiff(1:S, i);
                        Pi_i_d = zeros(r_dim, Mtot, nDeployments);
                        breve_g_i_d = zeros(r_dim, nDeployments);
                        for d_idx = 1:nDeployments
                            G_others = G_R(:, others, d_idx);      % Mtot x (S-1)
                            Pi_i_mat = null(G_others.').';          % r_dim x Mtot, orthonormal rows
                            Pi_i_d(:,:,d_idx) = Pi_i_mat;
                            breve_g_i_d(:,d_idx) = Pi_i_mat * G_R(:, i, d_idx);
                        end
                        Pi_all{i} = Pi_i_d;
                        breve_g_all{i} = breve_g_i_d;
                    end

                    for channel_db_idx = 1:size(channel_db_values,2)
                        channel_db_start = tic;
                        log_msg(verbosity_level, 2, 'Scheme = %s | Agent SNR = %d dB | Channel SNR = %g dB | S = %d', ...
                            scheme, agent_db, channel_db_values(scheme_idx,channel_db_idx,agent_db_idx), S);

                        if scheme == "EPC"
                            gamma_n = (dt/n_psd_constant) * Ps * E_mi_sqr * E_mag_g_sqr / db2magTen(channel_db_values(scheme_idx,channel_db_idx,agent_db_idx));
                        end

                        %% Set server noise power. Note the factor of 1/dt which is
                        % equivalent to applying an anti-aliasing filter.
                        scaled_n = sqrt( (gamma_n*n_psd_constant/dt) / 2 ) * n;
                        scaled_n_psd_constant = gamma_n * n_psd_constant; % == N0/2

                        %%% Expressions for original manuscript
                        V_arr_internal = @(w, ai, mi, ci) scaled_w_psd_constant .* sum(ci.^2 .* ai.^2 .* mi.^2, 2) .* mag_sqr_S0_internal(w) + scaled_n_psd_constant/2;
                        V_arr = @(w, ai, mi, ci) zeroIfnan(V_arr_internal(w, ai, mi, ci)) + (w == Bee | w == -Bee) .* (V_arr_internal(w-epsilon, ai, mi, ci) + V_arr_internal(w+epsilon, ai, mi, ci))./2;

                        H_arr_internal = @(w, ai, mi, ci) (1 + epsilon) ./ (V_arr_internal(w, ai, mi, ci) + epsilon);
                        H_arr = @(w, ai, mi, ci) zeroIfnan(H_arr_internal(w, ai, mi, ci)) + (w == Bee | w == -Bee) .* (H_arr_internal(w-epsilon, ai, mi, ci) + H_arr_internal(w+epsilon, ai, mi, ci))./2;

                        alpha_objective_array = @(ai, mi, ci) (2*pi) ./( ((sum(ci .* ai .* mi.^2)).^2) .* integral(@(w) (mag_sqr_S0(w)).^2 ./ V_arr(w, ai, mi, ci), -Inf, Inf, 'ArrayValued', true));
                        t0_objective_array = @(ai, mi, ci, arb_alpha) (2*pi) ./( (arb_alpha^2 .* (sum(ci .* ai .* mi.^2)).^2) .* integral(@(w) (w .* mag_sqr_S0(w)).^2 ./ V_arr(w,ai, mi, ci), -Inf, Inf, 'ArrayValued', true));

                        % --- Per-sensor known noise covariance eigenvalues for the residual test:
                        % K_i = m_i^2*(Nw/2)*D_template + (N0/2)/||breve_g_i||^2 * I (trial-independent).
                        Nw_over_2 = scaled_w_psd_constant;
                        N0_over_2 = scaled_n_psd_constant / (2*dt);
                        Ki_full_eig = cell(S,1);
                        for i = 1:S
                            a_i = reshape(mi_5d(1,i,1,1,:).^2 * Nw_over_2, 1, nDeployments);
                            b_i = N0_over_2 ./ sum(breve_g_all{i}.^2, 1);
                            Ki_full_eig{i} = lam_D .* a_i + b_i;   % K x nDeployments
                        end

                        y_no_attack = pagetranspose(out_cws_no_attack + scaled_n);   % SAME noise draw as every attacked y below

                        %% ---- Attack-independent quantities for this channel SNR (computed once) ----
                        if use_W
                            lambda_vals = (2/gamma_n) * lambda_B_base;               % 2M x nDeployments
                            W           = (1/sqrt(0.5*gamma_n)) * UB_transpose;      % 2M x 2M x nDeployments
                            W_bcast     = reshape(W, Mtot, Mtot, 1, nDeployments);
                            mu          = m.^2;
                            WG_Rmu      = reshape(pagemtimes(W, pagemtimes(G_R, reshape(mu, S, 1, nDeployments))), Mtot, nDeployments);

                            % Per-branch Wiener filters and their b_m^2-weighted sum
                            Qn_all = cell(Mtot, 1);
                            Q_sum  = zeros(K, K, nDeployments);
                            for m_idx = 1:Mtot
                                lambda_m = reshape(lambda_vals(m_idx,:), 1, 1, 1, nDeployments);
                                [~, Qn_all{m_idx}] = get_time_domain(Hm_arr(omega, lambda_m), dt, nDeployments);   % K x K x D
                                Q_sum = Q_sum + reshape(WG_Rmu(m_idx,:), 1, 1, nDeployments).^2 .* Qn_all{m_idx};
                            end

                            % Null geometries depend only on (nulled set, deployment) at this channel SNR.
                            geom_cache = containers.Map('KeyType','char','ValueType','any');

                            % --- Full-honest baseline: full S-sensor array, attacker never hijacked anyone
                            alpha_estimates_full_honest = zeros(1,1,nTrials,nDeployments);
                            t0_estimates_full_honest    = zeros(1,1,nTrials,nDeployments);
                            for d_idx = 1:nDeployments
                                geom = cached_null_geometry(geom_cache, [], mi_5d, g_tilde, gamma_w, gamma_n, Hm_arr, omega, dt, K, S, d_idx);
                                [alpha_b, t0_b, ~] = estimate_with_geometry(geom, y_no_attack, K, N, Tp, ...
                                    norm_fact, t, t0_true, mfTemplateFFT_raw, D_template, 1:nTrials, d_idx, dt, offset_idx, alpha_uses_true_t0);
                                alpha_estimates_full_honest(1,1,:,d_idx) = alpha_b;
                                t0_estimates_full_honest(1,1,:,d_idx) = t0_b;
                            end
                            full_honest_mse = alpha_t0_mse(alpha_estimates_full_honest, t0_estimates_full_honest, alpha_true, t0_true);

                            % --- CRLBs (2 x nDeployments: row 1 alpha, row 2 t0)
                            alpha_term = integral(@(w) sum(WG_Rmu.^2 .*          mag_sqr_S0(w).^2 ./ (lambda_vals.*mag_sqr_S0(w)+1), 1), -Inf, Inf, 'ArrayValued', true);
                            t0_term    = integral(@(w) sum(WG_Rmu.^2 .* (w.^2) .* mag_sqr_S0(w).^2 ./ (lambda_vals.*mag_sqr_S0(w)+1), 1), -Inf, Inf, 'ArrayValued', true);
                            crlb_pair = [2*pi ./ alpha_term; 2*pi ./ (alpha_true^2 .* t0_term)];

                            if attacker_enabled
                                lambda_vals_oracle = (2/gamma_n) * lambda_B_base_oracle;
                                W_oracle           = (1/sqrt(0.5*gamma_n)) * UB_transpose_oracle;
                                WG_Rmu_oracle      = reshape(pagemtimes(W_oracle, pagemtimes(G_R_oracle, reshape(mu_oracle, S_honest, 1, nDeployments))), Mtot, nDeployments);
                                alpha_term_oracle = integral(@(w) sum(WG_Rmu_oracle.^2 .*          mag_sqr_S0(w).^2 ./ (lambda_vals_oracle.*mag_sqr_S0(w)+1), 1), -Inf, Inf, 'ArrayValued', true);
                                t0_term_oracle    = integral(@(w) sum(WG_Rmu_oracle.^2 .* (w.^2) .* mag_sqr_S0(w).^2 ./ (lambda_vals_oracle.*mag_sqr_S0(w)+1), 1), -Inf, Inf, 'ArrayValued', true);
                                crlb_oracle_pair = [2*pi ./ alpha_term_oracle; 2*pi ./ (alpha_true^2 .* t0_term_oracle)];
                            end
                        else
                            ci = real(reshape(g_tilde, 1, S, 1, nDeployments));
                            crlb_pair = [reshape(alpha_objective_array(1, mi_4d, ci), 1, []);
                                         reshape(t0_objective_array(1, mi_4d, ci, alpha_true), 1, [])];
                        end

                        %% ---- Attack loop (innermost): only the compromised sensors' content changes ----
                        for pivot_idx = 1:numel(pivot_vals)
                            attack_type  = pivot_types(pivot_idx);   % attack mechanism
                            attack_label = pivot_vals(pivot_idx);    % label for logs and tables
                            attack_start = tic;

                            % Ground-truth compromised set for this run (empty when there is no attack).
                            if attack_type == "none"
                                true_atk = [];
                            else
                                true_atk = atk;
                            end

                            if attacker_enabled && attack_type ~= "none"
                                out_cws = out_cws_no_attack - out_cws_honest_atk + combine(ui_atk_by_attack{pivot_idx}, g_tilde(:,atk,:,:,:));
                            else
                                out_cws = out_cws_no_attack;
                            end
                            y = pagetranspose(out_cws + scaled_n);

                            %%% Baseline single-antenna case from previous manuscripts
                            if use_W == false
                                y_m = reshape(y, 1, K, nTrials, nDeployments);
                                mag_sqr_H = H_arr(omega, ai, mi_4d, ci);

                                % Get time domain matrix for Qn.
                                [~, Qn_matrix] = get_time_domain(mag_sqr_H,dt,nDeployments);

                                Omega_2_t0_search = reshape(sum(ui_t0_search .* ci, 2), K, 1, K, 1, nDeployments);
                                mf_with_y = dt*dt*pagemtimes(pagemtimes(reshape(real(y_m) + sum(ci .* bi, 2), 1, K, 1, nTrials, nDeployments), reshape(Qn_matrix,K,K,1,1,nDeployments)), Omega_2_t0_search);
                                [~,I] = max(mf_with_y, [], 3);

                                t0_estimates_for_plot = reshape((I-1)*dt, 1, 1, nTrials, nDeployments);

                                if alpha_uses_true_t0
                                    t0_estimates_for_alpha = t0_true*ones(1,1,nTrials,nDeployments);
                                else
                                    t0_estimates_for_alpha = t0_estimates_for_plot;
                                end

                                %%% Estimate alpha
                                Af = fft(mi_4d .* sensor_signal(t-t0_estimates_for_alpha, Tp, norm_fact), N, 1);
                                lag0 = round(Tp/dt);
                                Y_tensor =  dt*ifft(Af .* mfTemplateFFT, [], 1);
                                ui_tensor = Y_tensor(lag0+1 : (lag0+K), :, :, :);
                                Omega = sum(ui_tensor .* ci .* ai, 2);

                                resh_Qn = reshape(Qn_matrix,K,K,1,nDeployments);
                                resh_Omega = reshape(Omega,1,K,nTrials,nDeployments);

                                num = dt*dt*pagemtimes(pagemtimes(real(y_m) + sum(ci .* bi, 2),resh_Qn),pagetranspose(resh_Omega));
                                denom = dt*dt*pagemtimes(pagemtimes(resh_Omega,resh_Qn),pagetranspose(resh_Omega));

                                alpha_estimates = num ./ denom;

                            %%% Multi-antenna case
                            else
                                y_sq = reshape(y, K, num_antennas, nTrials, nDeployments);   % K x num_antennas x nTrials x nDeployments
                                y_R = permute(cat(2, real(y_sq), imag(y_sq)), [2 1 3 4]);   % 2M x K x nTrials x nDeployments

                                % --- Null-space isolation: recover each sensor's own contribution,
                                % including the attacker's, on the RAW received signal. ---
                                hat_u_all = zeros(K, S, nTrials, nDeployments);
                                for i = 1:S
                                    Pi_i_bcast = reshape(Pi_all{i}, r_dim, Mtot, 1, nDeployments);
                                    isolated = pagemtimes(Pi_i_bcast, y_R);   % r_dim x K x nTrials x nDeployments

                                    breve_g_i_bcast = reshape(breve_g_all{i}, 1, r_dim, 1, nDeployments);
                                    norm_sq_bcast = reshape(sum(breve_g_all{i}.^2, 1), 1, 1, 1, nDeployments);

                                    hat_u_i = pagemtimes(breve_g_i_bcast, isolated) ./ norm_sq_bcast;  % 1 x K x nTrials x nDeployments
                                    hat_u_all(:, i, :, :) = reshape(hat_u_i, K, 1, nTrials, nDeployments);
                                end

                                flagged = false(1, S, nTrials, nDeployments);

                                %% --- Per-sensor GLS delay and amplitude
                                t0_col_idx_per_sensor = zeros(S, nTrials, nDeployments);
                                t0_frac_per_sensor    = zeros(S, nTrials, nDeployments);   % sub-sample t0_i
                                c_scalar_at_best = zeros(S, nTrials, nDeployments);   % c_i(tau_hat)
                                G_scalar_at_best = zeros(S, nTrials, nDeployments);   % G_i(tau_hat)
                                Gd_scalar_at_best = zeros(S, nTrials, nDeployments);  % G'_i(tau_hat): derivative energy (t0_i precision)

                                for i = 1:S
                                    u_i_v_all = pagemtimes(U_D.', reshape(hat_u_all(:,i,:,:), K, nTrials, nDeployments));  % K x nTrials x nDeployments

                                    for d_idx = 1:nDeployments
                                        weights = 1 ./ Ki_full_eig{i}(:,d_idx);                 % K x 1

                                        % G_i(tau) for every candidate tau -- trial-independent
                                        G_i_all_tau = sum((D_template_v.^2) .* weights, 1);   % 1 x K

                                        % c_i(tau) for every candidate tau, all trials at once
                                        c_i_all_tau = (weights .* u_i_v_all(:,:,d_idx)).' * D_template_v;   % nTrials x K

                                        ratio = (c_i_all_tau.^2) ./ G_i_all_tau;                % nTrials x K

                                        [~, best_col] = max(ratio(:, offset_idx:end), [], 2);   % nTrials x 1
                                        best_col = best_col + offset_idx - 1;                   % column index into D_template

                                        lin_idx = sub2ind([nTrials, K], (1:nTrials).', best_col);
                                        c_scalar_at_best(i,:,d_idx) = c_i_all_tau(lin_idx);
                                        G_scalar_at_best(i,:,d_idx) = G_i_all_tau(best_col);
                                        Gd_i_all_tau = sum((D_deriv_v.^2) .* weights, 1);       % 1 x K
                                        Gd_scalar_at_best(i,:,d_idx) = Gd_i_all_tau(best_col);
                                        t0_col_idx_per_sensor(i,:,d_idx) = best_col;

                                        % Sub-sample t0_i: parabolic interpolation of the peak
                                        t0_frac_per_sensor(i,:,d_idx) = t(best_col) + refine_peak(ratio, best_col)*dt;
                                    end
                                end

                                t0_col_idx_per_trial = round(median(t0_col_idx_per_sensor, 1));   % 1 x nTrials x nDeployments

                                %% --- Normalized correlation vs. matched-filter template (median-delay) ---
                                corr_min = 0.95;
                                norm_corr_all = zeros(1, S, nTrials, nDeployments);
                                for d_idx = 1:nDeployments
                                    rho_ref_mat = D_template(:, reshape(t0_col_idx_per_trial(1,:,d_idx), 1, nTrials));   % K x nTrials
                                    rho_ref_norm = sqrt(sum(rho_ref_mat.^2, 1));            % 1 x nTrials
                                    for i = 1:S
                                        u_i = reshape(hat_u_all(:,i,:,d_idx), K, nTrials);   % K x nTrials
                                        num_r = sum(rho_ref_mat .* u_i, 1);
                                        den_r = sqrt(sum(u_i.^2, 1)) .* rho_ref_norm;
                                        norm_corr_all(1,i,:,d_idx) = reshape(num_r ./ den_r, 1, 1, nTrials);
                                    end
                                end

                                norm_corr_flagged = norm_corr_all < corr_min;
                                flagged = flagged | norm_corr_flagged;

                                %% --- Individual (per-sensor) t0 and alpha estimation ---
                                mi_sq = reshape(mi_5d.^2, S, 1, nDeployments);
                                t0_i_hat    = reshape(t0_frac_per_sensor, 1, S, nTrials, nDeployments);   % sub-sample
                                alpha_i_hat = reshape(c_scalar_at_best ./ (G_scalar_at_best .* mi_sq), 1, S, nTrials, nDeployments);

                                %% --- Robust outlier tests: flag if |x_i - median| / scale > k, where scale
                                % is a robust std estimate: 1.4826*MAD, Sn or Qn (all consistent for Gaussian data).
                                med_alpha = median(alpha_i_hat, 2);
                                mad_alpha = mad(alpha_i_hat, 1, 2);   % median absolute deviation across sensors
                                med_t0 = median(t0_i_hat, 2);
                                mad_t0 = mad(t0_i_hat, 1, 2);

                                mad_scale = 1.4826;   % MAD -> std for Gaussian data, so k is in "sigmas"
                                spread_names = ["mad", "sn", "qn"];
                                scale_alpha = {mad_scale*mad_alpha, robust_scale(alpha_i_hat, "sn"), robust_scale(alpha_i_hat, "qn")};
                                scale_t0    = {mad_scale*mad_t0,    robust_scale(t0_i_hat, "sn"),    robust_scale(t0_i_hat, "qn")};
                                z_alpha_all = cell(1, numel(spread_names));
                                z_t0_all    = cell(1, numel(spread_names));
                                for sp = 1:numel(spread_names)
                                    z_alpha_all{sp} = abs(alpha_i_hat - med_alpha) ./ scale_alpha{sp};
                                    z_t0_all{sp}    = abs(t0_i_hat    - med_t0)    ./ scale_t0{sp};
                                    z_alpha_all{sp}(isnan(z_alpha_all{sp})) = 0;   % 0/0: zero spread and the sensor sits at the median
                                    z_t0_all{sp}(isnan(z_t0_all{sp}))       = 0;
                                end
                                z_alpha = z_alpha_all{spread_names == detector_spread};   % scores used by the detector
                                z_t0    = z_t0_all{spread_names == detector_spread};

                                %%% With alpha
                                alpha_delta_k = 7;
                                alpha_flagged = z_alpha > alpha_delta_k;
                                flagged = flagged | alpha_flagged;

                                %%% With t0
                                t0_delta_k = 7;
                                t0_flagged = z_t0 > t0_delta_k;
                                flagged = flagged | t0_flagged;

                                %% --- ROC: count flags at every threshold on each test's own grid.
                                % Only the counts change with the thresholds -- nulling still uses the
                                % current corr_min / k above.
                                is_atk = false(1, S);
                                is_atk(true_atk) = true;
                                roc_idx = {pivot_idx, channel_db_idx, agent_db_idx, save_dim};
                                pick = @(x, mask) reshape(x(1, mask, :, :), [], 1);   % (sensor, trial, deployment) cases

                                % normcorr: flagged if norm_corr < corr_min
                                corr_atk = pick(norm_corr_all, is_atk);
                                corr_hon = pick(norm_corr_all, ~is_atk);
                                roc_corr_det(:, roc_idx{:}) = roc_corr_det(:, roc_idx{:}) + sum(corr_atk < roc_corr_grid, 1).';
                                roc_corr_fa(:, roc_idx{:})  = roc_corr_fa(:, roc_idx{:})  + sum(corr_hon < roc_corr_grid, 1).';

                                % alpha / t0 outlier tests, for every spread estimator: flagged if z > k
                                for sp = 1:numel(spread_names)
                                    roc_alpha_det(:, sp, roc_idx{:}) = roc_alpha_det(:, sp, roc_idx{:}) + sum(pick(z_alpha_all{sp}, is_atk)  > roc_k_grid, 1).';
                                    roc_alpha_fa(:, sp, roc_idx{:})  = roc_alpha_fa(:, sp, roc_idx{:})  + sum(pick(z_alpha_all{sp}, ~is_atk) > roc_k_grid, 1).';
                                    roc_t0_det(:, sp, roc_idx{:})    = roc_t0_det(:, sp, roc_idx{:})    + sum(pick(z_t0_all{sp}, is_atk)     > roc_k_grid, 1).';
                                    roc_t0_fa(:, sp, roc_idx{:})     = roc_t0_fa(:, sp, roc_idx{:})     + sum(pick(z_t0_all{sp}, ~is_atk)    > roc_k_grid, 1).';
                                end

                                roc_n_atk(roc_idx{:}) = roc_n_atk(roc_idx{:}) + numel(corr_atk);
                                roc_n_hon(roc_idx{:}) = roc_n_hon(roc_idx{:}) + numel(corr_hon);

                                %% --- Conventional baselines from the per-sensor estimates (independent of the detector)
                                % Rejection: conventional normalized-MAD rule, |x_i - median| / (1.4826*MAD) > conv_mad_k
                                z_mad_alpha = z_alpha_all{spread_names == "mad"};
                                z_mad_t0    = z_t0_all{spread_names == "mad"};
                                keep_a = z_mad_alpha <= conv_mad_k;   % 1 x S x nTrials x D
                                keep_t = z_mad_t0    <= conv_mad_k;
                                if conv_rejection == "joint"
                                    keep_a = keep_a & keep_t;
                                    keep_t = keep_a;
                                end
                                none_kept_a = ~any(keep_a, 2);        % 1 x 1 x nTrials x D -> fall back to the median
                                none_kept_t = ~any(keep_t, 2);
                                % Precision (1/variance) of each per-sensor estimate, for the weighted skipped mean:
                                %   alpha_i: G_i(tau_i) * m_i^4;   t0_i: G'_i(tau_i) * m_i^4 (the common alpha^2 cancels)
                                prec_alpha = reshape(G_scalar_at_best  .* mi_sq.^2, 1, S, nTrials, nDeployments);
                                prec_t0    = reshape(Gd_scalar_at_best .* mi_sq.^2, 1, S, nTrials, nDeployments);
                                conv_alpha = struct();
                                conv_t0    = struct();
                                for cm = conventional_methods
                                    switch cm
                                        case "median"
                                            a_est = med_alpha;
                                            t_est = med_t0;
                                        case "trimmed_mean"
                                            assert(2*conv_trim < S, 'conv_trim = %d trims every sensor (S = %d).', conv_trim, S);
                                            a_sorted = sort(alpha_i_hat, 2);
                                            t_sorted = sort(t0_i_hat, 2);
                                            a_est = mean(a_sorted(1, conv_trim+1:S-conv_trim, :, :), 2);
                                            t_est = mean(t_sorted(1, conv_trim+1:S-conv_trim, :, :), 2);
                                        case "skipped_mean"
                                            a_est = sum(alpha_i_hat .* keep_a, 2) ./ sum(keep_a, 2);
                                            t_est = sum(t0_i_hat    .* keep_t, 2) ./ sum(keep_t, 2);
                                        case "skipped_wmean"
                                            w_a = prec_alpha .* keep_a;
                                            w_t = prec_t0    .* keep_t;
                                            a_est = sum(w_a .* alpha_i_hat, 2) ./ sum(w_a, 2);
                                            t_est = sum(w_t .* t0_i_hat,    2) ./ sum(w_t, 2);
                                    end
                                    if startsWith(cm, "skipped")   % no sensor kept -> fall back to the median
                                        a_est(none_kept_a) = med_alpha(none_kept_a);
                                        t_est(none_kept_t) = med_t0(none_kept_t);
                                    end
                                    conv_alpha.(char(cm)) = a_est;   % 1 x 1 x nTrials x nDeployments
                                    conv_t0.(char(cm))    = t_est;
                                end

                                %% Nominal estimation using multiple receive antennas
                                % The t0 search and the alpha estimate are both linear in the per-branch
                                % filtered data, so build v = sum_m b_m*z_m*Qn_m once (Qn_m and Q_sum are
                                % precomputed above for this channel SNR).
                                z = pagemtimes(W_bcast, y_R);      % 2M x K x nTrials x nDeployments
                                v = zeros(nTrials, K, nDeployments);
                                for m_idx = 1:Mtot
                                    b_m = reshape(WG_Rmu(m_idx,:), 1, 1, nDeployments);
                                    v = v + b_m .* pagemtimes(permute(z(m_idx,:,:,:), [3 2 4 1]), Qn_all{m_idx});
                                end

                                mf_t0 = pagemtimes(v, D_template);
                                [~, I] = max(mf_t0(:, offset_idx:end, :), [], 2);
                                t0_col = I + offset_idx - 1;   % grid peak, nTrials x 1 x D
                                % Sub-sample t0_hat: parabolic interpolation of the peak (same step as the per-sensor t0_i)
                                t0_delta = zeros(size(t0_col));
                                for d_idx = 1:nDeployments
                                    t0_delta(:, 1, d_idx) = refine_peak(mf_t0(:, :, d_idx), t0_col(:, 1, d_idx));
                                end
                                t0_estimates_for_plot = reshape(t(t0_col) + t0_delta*dt, 1, 1, nTrials, nDeployments);

                                % Alpha given t0. Joint (default): rho(t - t0_hat) at the grid peak, a column of
                                % D_template (alpha uses the grid peak, as the per-sensor alpha_i does). Otherwise: the true t0.
                                if alpha_uses_true_t0
                                    Y_true = dt*ifft(fft(sensor_signal(t - t0_true, Tp, norm_fact), N, 1) .* mfTemplateFFT_raw);
                                    Rss    = repmat(Y_true(round(Tp/dt)+1 : round(Tp/dt)+K).', nTrials, 1, nDeployments);   % nTrials x K x D
                                else
                                    Rss    = permute(reshape(D_template(:, t0_col(:)), K, nTrials, nDeployments), [2 1 3]);
                                end
                                denom = dt*dt * sum(pagemtimes(Rss, Q_sum) .* Rss, 2);   % nTrials x 1 x D
                                alpha_estimates = reshape(dt*dt*sum(v .* Rss, 2) ./ denom, 1, 1, nTrials, nDeployments);

                                % Per-trial flagged sensor SET (face value -- no validation).
                                excluded_sets = cell(nTrials, nDeployments);
                                for d_idx = 1:nDeployments
                                    for tr = 1:nTrials
                                        excluded_sets{tr,d_idx} = find(flagged(1,:,tr,d_idx));   % row vec, possibly empty
                                    end
                                end

                                if attacker_enabled
                                    alpha_estimates_undefended = alpha_estimates;
                                    t0_estimates_undefended    = t0_estimates_for_plot;

                                    % --- Oracle ceiling: null the TRUE attacker_idx directly
                                    alpha_estimates_oracle = zeros(size(alpha_estimates));
                                    t0_estimates_oracle    = zeros(size(t0_estimates_for_plot));
                                    for d_idx = 1:nDeployments
                                        geom = cached_null_geometry(geom_cache, true_atk, mi_5d, g_tilde, gamma_w, gamma_n, Hm_arr, omega, dt, K, S, d_idx);
                                        [alpha_o, t0_o, ~] = estimate_with_geometry(geom, y, K, N, Tp, ...
                                            norm_fact, t, t0_true, mfTemplateFFT_raw, D_template, 1:nTrials, d_idx, dt, offset_idx, alpha_uses_true_t0);
                                        alpha_estimates_oracle(1,1,:,d_idx) = alpha_o;
                                        t0_estimates_oracle(1,1,:,d_idx) = t0_o;
                                    end
                                end

                                alpha_final = alpha_estimates;
                                t0_final = t0_estimates_for_plot;

                                %% --- Face value validation: group trials by exact flagged set, null once per group
                                tp_total = 0;
                                fp_total = 0;
                                fn_total = 0;

                                for d_idx = 1:nDeployments
                                    [unique_sets, group_idx] = unique_cell_sets(excluded_sets(:,d_idx));

                                    for grp = 1:numel(unique_sets)
                                        A = unique_sets{grp};
                                        trial_members = find(group_idx == grp);
                                        n_members = numel(trial_members);

                                        % Count detection outcomes for every group, including "nothing flagged".
                                        tp_total = tp_total + numel(intersect(A, true_atk)) * n_members;
                                        fp_total = fp_total + numel(setdiff(A, true_atk))   * n_members;
                                        fn_total = fn_total + numel(setdiff(true_atk, A))   * n_members;

                                        if isempty(A)
                                            continue   % nothing flagged -- keep undefended estimate
                                        end

                                        geom = cached_null_geometry(geom_cache, A, mi_5d, g_tilde, gamma_w, gamma_n, Hm_arr, omega, dt, K, S, d_idx);
                                        [alpha_A, t0_A, ~] = estimate_with_geometry(geom, y, K, N, Tp, ...
                                            norm_fact, t, t0_true, mfTemplateFFT_raw, D_template, trial_members, d_idx, dt, offset_idx, alpha_uses_true_t0);

                                        alpha_final(1,1,trial_members,d_idx) = alpha_A;
                                        t0_final(1,1,trial_members,d_idx) = t0_A;
                                    end
                                end

                                log_msg(verbosity_level, 4, '[%s] Detection accuracy: TP=%d, FP=%d, FN=%d (%d null geometries cached)', ...
                                    attack_label, tp_total, fp_total, fn_total, geom_cache.Count);

                                % --- Honest-sensor false-flag rate per test
                                honest_idx_diag = setdiff(1:S, true_atk);
                                fp_stats = containers.Map('KeyType','char','ValueType','double');
                                fp_stats('alpha_i')  = mean(alpha_flagged(1,honest_idx_diag,:,:), 'all');
                                fp_stats('t0_i')     = mean(t0_flagged(1,honest_idx_diag,:,:), 'all');
                                fp_stats('normcorr') = mean(norm_corr_flagged(1,honest_idx_diag,:,:), 'all');
                                fp_stats('combined') = mean(flagged(1,honest_idx_diag,:,:), 'all');
                                fp_report(fp_stats, verbosity_level);

                                det_idx = {pivot_idx, channel_db_idx, agent_db_idx, save_dim};
                                det_tp(det_idx{:}) = det_tp(det_idx{:}) + tp_total;   % summed over batches
                                det_fp(det_idx{:}) = det_fp(det_idx{:}) + fp_total;
                                det_fn(det_idx{:}) = det_fn(det_idx{:}) + fn_total;
                                fp_rate_combined(det_idx{:}) = fp_rate_combined(det_idx{:}) + fp_stats('combined') * nDeployments / nDeployments_total;

                                alpha_estimates = alpha_final;
                                t0_estimates_for_plot = t0_final;
                            end

                            %% Store this attack's results
                            slot = {1, agent_db_idx, channel_db_idx, ':', scheme_idx, pivot_idx, save_dim, dep_range};
                            rho_empirical_var(slot{:}) = [reshape(var(alpha_estimates,0,3),1,[]); reshape(var(t0_estimates_for_plot,0,3),1,[])];
                            rho_empirical_mse(slot{:}) = alpha_t0_mse(alpha_estimates, t0_estimates_for_plot, alpha_true, t0_true);
                            rho_crlb(slot{:}) = crlb_pair;

                            if attacker_enabled && use_W
                                rho_empirical_mse_full_honest(slot{:}) = full_honest_mse;
                                rho_empirical_mse_undefended(slot{:})  = alpha_t0_mse(alpha_estimates_undefended, t0_estimates_undefended, alpha_true, t0_true);
                                rho_empirical_mse_oracle(slot{:})      = alpha_t0_mse(alpha_estimates_oracle,     t0_estimates_oracle,     alpha_true, t0_true);
                                for cm = conventional_methods
                                    rho_mse_conv.(char(cm))(slot{:}) = alpha_t0_mse(conv_alpha.(char(cm)), conv_t0.(char(cm)), alpha_true, t0_true);
                                end
                                if attack_type == "none"
                                    rho_crlb_oracle(slot{:}) = crlb_pair;          % no one to remove: oracle bound = full-array bound
                                else
                                    rho_crlb_oracle(slot{:}) = crlb_oracle_pair;
                                end
                            end

                            attack_time = toc(attack_start);
                            log_msg(verbosity_level, 3, 'Attack "%s": %dm %.1fs', attack_label, floor(attack_time/60), mod(attack_time,60));
                        end

                        % Checkpoint after each channel-SNR point so a crash doesn't lose finished results.
                        save(experiment + "_checkpoint.mat", "rho_*", "det_*", "roc_*", "fp_rate_combined", "pivot_vals", ...
                            "batch_idx", "agent_db_idx", "sensor_idx", "scheme_idx", "channel_db_idx");

                        total_channel_db_time = toc(channel_db_start);
                        log_msg(verbosity_level, 3, 'Runtime for channel SNR = %g dB (all attacks): %dm %.1fs\n', ...
                            channel_db_values(scheme_idx,channel_db_idx,agent_db_idx), floor(total_channel_db_time/60), mod(total_channel_db_time,60));
                    end
                end
            end
        end

        % Batch runtime, with an estimate of the time left from the average batch so far
        batch_time(batch_idx) = toc(batch_start);
        eta = mean(batch_time(1:batch_idx)) * (nBatches - batch_idx);
        log_msg(verbosity_level, 1, 'Batch %d of %d done in %dm %.1fs (%.2f s per deployment) | elapsed %dm %.1fs | est. remaining %dm %.1fs', ...
            batch_idx, nBatches, floor(batch_time(batch_idx)/60), mod(batch_time(batch_idx),60), batch_time(batch_idx)/nDeployments, ...
            floor(toc(loopTic)/60), mod(toc(loopTic),60), floor(eta/60), mod(eta,60));
    end   % deployment batches

    % Stop run timer.
    runtime = toc(loopTic);
    % Display total run time.
    log_msg(verbosity_level, 1, '%s took %dm %.1fs', experiment, floor(runtime/60), mod(runtime,60));
    
    %% Compute averages across deployments.
    if nDeployments_total == 1
        avg_dep_var = rho_empirical_var;
        avg_dep_mse = rho_empirical_mse;
        avg_dep_crlb = rho_crlb;
        avg_dep_bias = rho_empirical_bias;
        avg_dep_mse_undefended = rho_empirical_mse_undefended;
        avg_dep_mse_oracle = rho_empirical_mse_oracle;   % (or mean(...) branch, matching the existing if/else)
        avg_dep_mse_full_honest = rho_empirical_mse_full_honest;
        % avg_dep_mse_basis_supp = rho_empirical_mse_basis_supp;
        avg_dep_crlb_oracle = rho_crlb_oracle;
    else
        avg_dep_var = mean(rho_empirical_var,length(size(rho_empirical_var)));
        avg_dep_mse = mean(rho_empirical_mse,length(size(rho_empirical_mse)));
        avg_dep_crlb = mean(rho_crlb,length(size(rho_crlb)));
        avg_dep_bias = mean(rho_empirical_bias,length(size(rho_empirical_bias)));
        avg_dep_mse_undefended = mean(rho_empirical_mse_undefended,length(size(rho_empirical_mse_undefended)));
        avg_dep_mse_oracle = mean(rho_empirical_mse_oracle,length(size(rho_empirical_mse_oracle)));   % (or mean(...) branch, matching the existing if/else)
        avg_dep_mse_full_honest = mean(rho_empirical_mse_full_honest,length(size(rho_empirical_mse_full_honest)));
        % avg_dep_mse_basis_supp = mean(rho_empirical_mse_basis_supp,length(size(rho_empirical_mse_basis_supp)));
        avg_dep_crlb_oracle = mean(rho_crlb_oracle,length(size(rho_crlb_oracle)));
    end

    avg_dep_mse_conv = struct();   % conventional baselines, one field per method
    for cm = conventional_methods
        avg_dep_mse_conv.(char(cm)) = mean(rho_mse_conv.(char(cm)), 8);   % dim 8 = deployments
    end

    txt_list = ["undefended", "oracle", "full_honest"]; % , "basis_supp"];
    save_images = false;

    files_to_save = ["avg_dep_mse_" + txt_list, "avg_dep_var", "avg_dep_mse", "avg_dep_bias", "avg_dep_crlb", "avg_dep_crlb_oracle",...
        "agent_db_values", "channel_snr", "selected_schemes", "sensor_vals", "channel_db_values",...
        "experiment", "sensor_dimension", "save_images", "K", "nDeployments", "nDeployments_total", "constraints", "transforms", "coeff_types", "approaches",...
        "iter_arr", "iter_length", "all_strats", "num_strats", "save_images", "alpha_uses_true_t0",...
        "T0", "Tp", "alpha_true", "t0_true", "pivot_vals", "Mtot", "N", "S", "S_max","attacks", "attacker_enabled", "pivot_types", "pivot_shift", "batch_time", "runtime", ...
        "avg_dep_mse_conv", "conventional_methods", "conv_trim", "conv_mad_k", "conv_rejection"];
    
    save(experiment + "_results.mat", files_to_save{:})
    if experiment == "time_shift_sensitivity"
        save(experiment + "_results.mat", "ts_delta", "ts_multiples", "ts_shift_s", "ts_shift_pct", "-append");
    end

    %% Compile per-attack summary (first strategy / scheme / agent SNR / sensor count)
    method_names = ["defended", "undefended", "oracle", "full_honest", conventional_methods];
    method_data  = [{avg_dep_mse, avg_dep_mse_undefended, avg_dep_mse_oracle, avg_dep_mse_full_honest}, ...
        cellfun(@(f) avg_dep_mse_conv.(f), cellstr(conventional_methods), 'UniformOutput', false)];
    rows = {};
    for p = 1:numel(pivot_vals)
        for c = 1:size(channel_db_values, 2)
            row = {pivot_vals(p), channel_db_values(1,c,1)};
            for mth = 1:numel(method_data)
                row = [row, {method_data{mth}(1,1,c,1,1,p,1), method_data{mth}(1,1,c,2,1,p,1)}];
            end
            rows(end+1,:) = [row, {det_tp(p,c,1,1), det_fp(p,c,1,1), det_fn(p,c,1,1), fp_rate_combined(p,c,1,1)}];
        end
    end
    var_names = ["attack", "channel_snr_db", reshape(method_names + ["_mse_alpha"; "_mse_t0"], 1, []), "TP", "FP", "FN", "honest_fp_rate"];
    summary_tbl = cell2table(rows, 'VariableNames', cellstr(var_names));
    disp(summary_tbl)
    writetable(summary_tbl, experiment + "_summary.csv");
    save(experiment + "_results.mat", "summary_tbl", "det_tp", "det_fp", "det_fn", "fp_rate_combined", ...
        "roc_corr_grid", "roc_k_grid", "roc_spread_names", "detector_spread", "roc_corr_det", "roc_corr_fa", ...
        "roc_alpha_det", "roc_alpha_fa", "roc_t0_det", "roc_t0_fa", "roc_n_atk", "roc_n_hon", ...
        "corr_min", "alpha_delta_k", "t0_delta_k", "mad_scale", "-append");   % current operating point
    
    %% Plotting placeholder
end
