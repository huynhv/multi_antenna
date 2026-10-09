%%% 10-6-2026 --> revisit tomorrow

% Clear all variables and close all existing figures.
clearvars
close all

addpath('.\functions')
verbosity_level = 5;   % 1 = major sections only, 2 = + sub-steps/channel-SNR, 3 = + diagnostics

% Set seed for reproducibility
seed = 22;
rng(seed);

use_W = true;

nDeployments = 10;
% Set the number of measurements per sensor, starting from 0
K = 301;
% Set number of trials.
nTrials = 500;

% Set number of sensors
sensor_vals = 9; % [5,10]
dropout_vals = 0;
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
t0_true = 1.6;

% Estimation mode: false = joint (alpha evaluated at the estimated t0),
%                  true  = alpha evaluated at the true t0.
alpha_uses_true_t0 = false;

% Set the number of antennas
num_antennas = 9;
Mtot = 2*num_antennas;

%%% Always generate w, n, and channel gain for maximum S, then use subsets for different S trials
S_max = max(sensor_vals);
% Generate sensor noise with dimensions: K x S x trials x deployments.
all_w = randn(K,S_max,1,nTrials,nDeployments);
% Generate server noise.
n = (randn(K,1,num_antennas,nTrials,nDeployments) + 1i*randn(K,1,num_antennas,nTrials,nDeployments));
% Generate channel gains.
all_gi = (randn(1,S_max,num_antennas,1,nDeployments) + 1i*randn(1,S_max,num_antennas,1,nDeployments))/sqrt(2);

%% Generate mi and ti
max_ti = 1;
min_ti = 0.8; 
path_loss_exp = 0.5;
shadow_sigma_dB = 0.2;

m_ref = 1;              % reference amplitude at distance min_ti
c_gain = m_ref * min_ti^path_loss_exp;   % anchors: unshadowed mi = m_ref at d = min_ti

all_di = sort(unifrnd(min_ti, max_ti, 1, S_max, 1, 1, nDeployments), 2, 'ascend');
all_ti = all_di;   % v = 1, per earlier discussion

shadowing = 10.^(shadow_sigma_dB/20 * randn(1,S_max,1,1,nDeployments));
all_mi = c_gain ./ (all_di.^path_loss_exp) .* shadowing;

%%
selected_schemes = "EPC";

% Define Rayleigh distribution parameters.
rayleigh_factor = 1/sqrt(2);
E_mag_g_sqr = 2*rayleigh_factor^2;
E_mi_sqr = mean(all_mi(:).^2); % --> could generate empirical value here

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
attacker_idx     = [4,5,6]; % [2,4]      % sensor index/indices (within 1:S) that are compromised
attacker_db      = 0;      % attacker "SNR", same convention as agent_db (see below)
attacker_seed = 42;

w_psd_constant = 1;
n_psd_constant = 1;
N = 2*K - 1;

%% Subspace-projection residual test: honest-signal basis Phi and its
t0_grid_top = repmat(reshape(t,1,1,[]), 1,1,1,nDeployments);
[~, R00_full] = mf_integral_fft(sensor_signal(t-t0_grid_top,Tp,norm_fact), sensor_signal(t,Tp,norm_fact), 1, 1, K, dt, Tp);
D_template = reshape(R00_full(:,:,:,1), K, K);   % K x K, [D_template]_{k,l} = rho(t_k - t_l)

[U_D, Sigma_D, ~] = svd(D_template, 'econ');
sv_energy = cumsum(diag(Sigma_D).^2) / sum(diag(Sigma_D).^2);
d_sub = find(sv_energy >= 0.999, 1, 'first');
Phi = U_D(:, 1:d_sub);           % K x d_sub: honest-signal subspace

% Eigenvalues of D_template itself (already computed via U_D/Sigma_D above)
lam_D = diag(Sigma_D);   % K x 1
% Project D_template's columns into U_D's eigenbasis once -- reused every
% (sensor, deployment) below since it doesn't depend on either.
D_template_v = U_D.' * D_template;   % K x K
mfTemplateFFT_raw = fft(sensor_signal(t, Tp, norm_fact), N, 1);   % used by estimate_with_geometry

log_msg(verbosity_level, 1, 'Subspace-projection residual test: d = %d, K-d = %d', d_sub, K-d_sub);

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
experiment_list = "multi_antenna_disjoint";

for experiment_idx = 1:numel(experiment_list)
    experiment = experiment_list(experiment_idx);

    pad_str = '----------';
    msg = [pad_str  ' ' 'Running: ' char(experiment) ' ' pad_str];
    log_msg(verbosity_level, 1, '%s', msg);

    % Configure strategies
    agent_db_values = 0;
    channel_snr = linspace(0,10,5);

    % "noise" | "sawtooth" | "sinusoid" | "peak" | "flip" | "time shift" | "amplitude scale"
    attacks = ["noise"];

    if numel(attacks) == 0
        error("No attack type has been specified!")
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
    rho_crlb = zeros(num_strats,length(agent_db_values),size(channel_db_values,2),2,length(selected_schemes),numel(pivot_vals),sensor_dimension,nDeployments);
    rho_crlb_oracle = rho_crlb;   % Theoretical CRLB with the attacker's channel treated as absent -- upper bound, only under attack
    
    rho_empirical_bias = rho_crlb;
    rho_empirical_var = rho_crlb;
    rho_empirical_mse = rho_crlb;
    rho_empirical_mse_baseline = rho_crlb;   % Full S-sensor array, attacker never hijacked anyone -- "what if there had been no attack at all"
    rho_empirical_mse_undefended = rho_crlb;   % MSE with NO defense applied -- only meaningful/populated under attack
    rho_empirical_mse_oracle = rho_crlb;   % MSE with attacker PERFECTLY, freely removed -- upper bound, only under attack
    rho_empirical_mse_bs = rho_crlb;
    rho_empirical_mse_trimmed = rho_crlb;

    % Start run timer
    log_msg(verbosity_level, 1, 'Starting runtime...');
    loopTic = tic;
    
    % Iterate through agent snr values.
    for pivot_idx = 1:length(pivot_vals)
        attack_type = pivot_vals(pivot_idx);

        rng(attacker_seed);
        f_nominal = Bee/(2*pi);   % = 1/(2*Tp), reference frequency near the honest passband

        switch attack_type
            case "noise"
                all_attacker_noise = randn(K,S_max,1,nTrials,nDeployments);

            case "sinusoid"
                f_c   = f_nominal * (0.5 + rand(1,S_max));
                phase = 2*pi*rand(1,S_max);
                waveform = sin(2*pi*t.*f_c + phase);               % K x S_max
                all_attacker_noise = repmat(reshape(waveform,K,S_max,1,1,1), 1,1,1,nTrials,nDeployments);

            case "peak"
                peak_same_location = true;
                candidate_locs = [dt, Tp];   % the two worst-case candidates identified earlier

                if peak_same_location
                    loc = candidate_locs(randi(2));
                    peak_loc = repmat(loc, 1, S_max);
                else
                    peak_loc = candidate_locs(randi(2, 1, S_max));
                end

                waveform = zeros(K, S_max);
                for a = 1:S_max
                    [~, idx] = min(abs(t - peak_loc(a)));
                    waveform(idx, a) = 1;
                end
                all_attacker_noise = repmat(reshape(waveform,K,S_max,1,1,1), 1,1,1,nTrials,nDeployments);

            case "flip"
                all_attacker_noise = [];   % unused placeholder

            case "amplitude scale"
                all_attacker_noise = [];   % unused placeholder

            case "time shift"
                all_attacker_noise = [];   % unused placeholder

            otherwise
                error("Unknown attack_type: %s", attack_type);
        end

        %% Normalize EVERY attack type to unit total (dt-weighted) energy per slot,
        if ~isempty(all_attacker_noise)
            waveform_energy = squeeze(sum(all_attacker_noise(:,:,1,1,1).^2, 1)) * dt;   % 1 x S_max
            norm_factor = reshape(1./sqrt(waveform_energy), 1, S_max, 1, 1, 1);
            all_attacker_noise = all_attacker_noise .* norm_factor;
        end

        rng(seed);

        for agent_db_idx = 1:length(agent_db_values)
            agent_db = agent_db_values(agent_db_idx);
        
            % Set sensor noise power. Note the factor of 1/dt, which is equivalent to
            % applying an anti-aliasing filter.
            gamma_w = (dt/w_psd_constant) * (Ps / db2magTen(agent_db));
            scaled_w = sqrt(gamma_w * w_psd_constant /dt) * all_w; % --> variance of this should be gamma_w/dt
            scaled_w_psd_constant = gamma_w * w_psd_constant; %% = Nw/2
    
            % Attacker noise power, using the SAME convention as gamma_w above --
            % just substituting attacker_db for agent_db.
            gamma_w_attacker = (dt/w_psd_constant) * (Ps / db2magTen(attacker_db));
            scaled_attacker_noise = sqrt(gamma_w_attacker * w_psd_constant / dt) * all_attacker_noise;

            % Iterate through sensor values.
            for sensor_idx = 1:length(sensor_vals)
                S = sensor_vals(sensor_idx);
        
                % Iterate through schemes.
                for scheme_idx = 1:length(selected_schemes)
                    scheme = selected_schemes(scheme_idx);
        
                    % Select subset of mi, ti, and gi values.
                    mi_5d = all_mi(:,1:S,:,:,:);
                    mi_4d = reshape(mi_5d, 1, S, 1, nDeployments);
                    ti = all_ti(:,1:S,:,:,:);
                    gi = all_gi(:,1:S,:,:,:);

                    noisy_unified_xi = alpha_true .* (mi_5d .* sensor_signal(t-t0_true, Tp, norm_fact)) + scaled_w(:,1:S,:,:,:);
                    % Compute ui here so we don't have to recompute FFT for each strat
                    [~, ui] = mf_integral_fft(noisy_unified_xi, mi_5d .* sensor_signal(t, Tp, norm_fact), 1, 1, K, dt, Tp);

                    % Capture the PRISTINE, never-attacked ui
                    ui_no_attack = ui;

                    %% --- Byzantine attacker: replace compromised sensor(s)' transmitted ui ---
                    if attacker_enabled
                        atk = attacker_idx(attacker_idx <= S);
                        switch attack_type
                            case "flip"              % -alpha*m_i^2*rho(t-t0) - w_i(t)
                                ui(:,atk,:,:,:) = -ui(:,atk,:,:,:);
                            case "amplitude scale"   % beta*(alpha*m_i^2*rho(t-t0) + w_i(t))
                                attacker_beta = 0.1;
                                ui(:,atk,:,:,:) = attacker_beta * ui(:,atk,:,:,:);
                            case "time shift"        % shifted copy, truncated to [0, T0]
                                t0_attacker = 3;
                                assert(abs((t0_attacker - t0_true)/dt - round((t0_attacker - t0_true)/dt)) < 1e-9, 'time shift is not a multiple of dt');
                                ui(:,atk,:,:,:) = shift_waveform(ui(:,atk,:,:,:), t, t0_attacker - t0_true);
                            otherwise
                                ui(:,atk,:,:,:) = scaled_attacker_noise(:,atk,:,:,:);
                        end
                    end

                    %% Define search space for t0
                    if ~use_W
                        t0_search_space = repmat(reshape(t,1,1,[]), 1,1,1,nDeployments);
                        [~, ui_t0_search] = mf_integral_fft(mi_4d .* sensor_signal(t-t0_search_space, Tp, norm_fact), mi_4d .* sensor_signal(t, Tp, norm_fact), 1, 1, K, dt, Tp);
                        mfTemplateFFT = fft(mi_4d .* sensor_signal(t, Tp, norm_fact), N, 1);
                    end

                    if scheme == "EPC"
                        % Apply exact phase compensation
                        ref_ant_idx = 1;
                        g_ref = gi(:, :, ref_ant_idx, :, :);
                        phase_comp = exp(-1j * angle(g_ref));
                        g_tilde = gi .* phase_comp;
                    end
                        
                    out_cws = sum(g_tilde .* ui, 2);
                    out_cws_no_attack = sum(g_tilde .* ui_no_attack, 2);

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
                    % reused across every channel-SNR point and trial for this S/scheme) ---
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

                    % Iterate through dropout values.
                    for dropout_idx = 1:length(dropout_vals)
                        save_dim = sensor_idx;

                        % Iterate through channel snr values.
                        for channel_db_idx = 1:size(channel_db_values,2)
                            channel_db_start = tic;
    
                            % Print statement for at-a-glance performance.
                            log_msg(verbosity_level, 2, 'Scheme %s | Attack Type = %s | Agent SNR = %d dB | Channel SNR = %g dB | S = %d | Dropout = %d', ...
                                scheme, attack_type, agent_db_values(agent_db_idx), channel_db_values(scheme_idx,channel_db_idx,agent_db_idx), S, dropout_vals(dropout_idx));

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

                            % --- Per-sensor known noise covariance for the residual test:
                            % K_i = m_i^2*(Nw/2)*Rss_Psi + (N0/2)/||breve_g_i||^2 * I.
                            % Computed once per channel-SNR point (deterministic given known
                            % parameters -- does NOT depend on trial noise realizations, so
                            % it's reused across all nTrials below). ---
                            Nw_over_2 = scaled_w_psd_constant;   % == Nw/2
                            N0_over_2 = scaled_n_psd_constant / (2*dt);   % == N0/2
                            
                            % Ki_diag_all = cell(S,1);
                            Ki_full_eig = cell(S,1);
                            for i = 1:S
                                a_i = reshape(mi_5d(1,i,1,1,:).^2 * Nw_over_2, 1, nDeployments);
                                norm_sq_i = zeros(1,nDeployments);
                                for d_idx = 1:nDeployments
                                    norm_sq_i(d_idx) = sum(breve_g_all{i}(:,d_idx).^2);
                                end
                                b_i = N0_over_2 ./ norm_sq_i;

                                % Ki_diag_all{i} = 1 ./ (lam_Rss .* a_i + b_i);   % (K-d_sub) x nDeployments
                                Ki_full_eig{i} = lam_D .* a_i + b_i;            % K x nDeployments (eigenVALUES, not inverted -- see usage below)
                            end

                            %%%%%%%%%%%%%%%%% ESTIMATION %%%%%%%%%%%%%%%%%
                            total_est_start = tic;
    
                            log_msg(verbosity_level, 3, 'Begin ML estimation');
                            for strat_idx = 1:num_strats
                                strat = all_strats(strat_idx);
                                log_msg(verbosity_level, 4, 'Running ML estimation for %s', strat);
        
                                ai = 1;
                                bi = 0;

                                y = pagetranspose(out_cws + scaled_n);
                                y_no_attack = pagetranspose(out_cws_no_attack + scaled_n);   % SAME noise draw, only the attacked sensor's content differs

                                %%% Baseline single-antenna case from previous manuscripts
                                if use_W == false
                                    ci = real(reshape(g_tilde, 1, S, 1, nDeployments));
                                    y_m = reshape(y, 1, K, nTrials, nDeployments);
                                    mag_sqr_H = H_arr(omega, ai, mi_4d, ci);
    
                                    % Get time domain matrix for Qn.
                                    [~, Qn_matrix] = get_time_domain(mag_sqr_H,dt,nDeployments);

                                    Omega_2_t0_search = reshape(sum(ui_t0_search .* ci .* ai, 2), K, 1, K, 1, nDeployments);
                                    mf_with_y = dt*dt*pagemtimes(pagemtimes(reshape(real(y_m) + sum(ci .* bi, 2), 1, K, 1, nTrials, nDeployments), reshape(Qn_matrix,K,K,1,1,nDeployments)), Omega_2_t0_search);
                                    [max_y_vals,I] = max(mf_with_y, [], 3);

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
                                    % ---- channel-DEPENDENT: rescale W and lambda (per channel SNR) ----
                                    lambda_vals = (2/gamma_n)      * lambda_B_base;   % 2M x nDeployments
                                    W           = (1/sqrt(0.5*gamma_n)) * UB_transpose;  % 2M x 2M x nDeployments

                                    if attacker_enabled
                                        lambda_vals_oracle = (2/gamma_n) * lambda_B_base_oracle;
                                        W_oracle           = (1/sqrt(0.5*gamma_n)) * UB_transpose_oracle;
                                        WG_Rmu_oracle      = reshape(pagemtimes(W_oracle, pagemtimes(G_R_oracle, reshape(mu_oracle, S_honest, 1, nDeployments))), Mtot, nDeployments);
                                    end

                                    y_sq = reshape(y, K, num_antennas, nTrials, nDeployments);   % K x num_antennas x nTrials x nDeployments
                                    % g_sq = reshape(g_tilde, S, num_antennas, nDeployments);   % S x num_antennas x nDeployments

                                    y_R = cat(2, real(y_sq), imag(y_sq));   % K x 2M x nTrials x nDeployments
                                    y_R = permute(y_R, [2 1 3 4]);   % 2M x K x nTrials x nDeployments

                                    y_R_bs = pagemtimes(pagemtimes(y_R, Phi), Phi.');   % project onto honest subspace, same K dims

                                    % --- Null-space isolation: recover each sensor's own
                                    % contribution, including the attacker's, on the RAW
                                    % (unprojected, un-decorrelated) received signal. ---
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
                                    % honest_idx = setdiff(1:S, attacker_idx);
                                    
                                    %% --- Practical hat_t0 for the sign-flip check: robust
                                    % version. (1) EXCLUDES sensors already flagged by T_i
                                    % corr_per_sensor = dt * pagemtimes(D_template.', reshape(hat_u_all, K, S*nTrials*nDeployments));   % K x (S*nTrials*nDeployments)
                                    % [~, best_col_per_sensor] = max(corr_per_sensor(offset_idx:end,:), [], 1);
                                    % 
                                    % best_col_per_sensor = best_col_per_sensor + offset_idx - 1;
                                    % t0_col_idx_per_sensor = reshape(best_col_per_sensor, S, nTrials, nDeployments) - 1;   % per-sensor delay guess

                                    t0_col_idx_per_sensor = zeros(S, nTrials, nDeployments);
                                    c_scalar_at_best = zeros(S, nTrials, nDeployments);   % cache c_i(tau_hat) for reuse below
                                    G_scalar_at_best = zeros(S, nTrials, nDeployments);   % cache G_i(tau_hat) for reuse below

                                    for i = 1:S
                                        u_i_v_all = pagemtimes(U_D.', reshape(hat_u_all(:,i,:,:), K, nTrials, nDeployments));  % K x nTrials x nDeployments

                                        for d_idx = 1:nDeployments
                                            eig_i_d = Ki_full_eig{i}(:,d_idx);        % K x 1
                                            weights = 1 ./ eig_i_d;                    % K x 1

                                            % G_i(tau) for every candidate tau -- trial-independent, computed once
                                            G_i_all_tau = sum((D_template_v.^2) .* weights, 1);   % 1 x K

                                            % c_i(tau) for every candidate tau, all trials at once
                                            u_i_v_d = u_i_v_all(:,:,d_idx);                        % K x nTrials
                                            c_i_all_tau = (weights .* u_i_v_d).' * D_template_v;   % nTrials x K

                                            ratio = (c_i_all_tau.^2) ./ G_i_all_tau;                % nTrials x K, broadcast over rows

                                            [~, best_col] = max(ratio(:, offset_idx:end), [], 2);   % nTrials x 1
                                            best_col = best_col + offset_idx - 1;                   % 1-based column index into D_template

                                            lin_idx = sub2ind([nTrials, K], (1:nTrials).', best_col);
                                            c_scalar_at_best(i,:,d_idx) = c_i_all_tau(lin_idx);
                                            G_scalar_at_best(i,:,d_idx) = G_i_all_tau(best_col);

                                            t0_col_idx_per_sensor(i,:,d_idx) = best_col;
                                        end
                                    end

                                    t0_col_idx_per_trial = round(median(t0_col_idx_per_sensor, 1));   % 1 x nTrials x nDeployments

                                    %% --- Normalized correlation vs. matched-filter template (median-delay) ---
                                    % r_i = <u_i, rho_ref> / (||u_i|| ||rho_ref||), cosine similarity in [-1,1].
                                    % rho_ref = D_template(:, t0_col_idx_per_trial), the median-delay template.
                                    corr_min = 0.95;
                                    norm_corr_all = zeros(1, S, nTrials, nDeployments);
                                    for d_idx = 1:nDeployments
                                        % K x nTrials matrix of templates, one column per trial (median delay)
                                        rho_ref_mat = D_template(:, reshape(t0_col_idx_per_trial(1,:,d_idx), 1, nTrials));
                                        rho_ref_norm = sqrt(sum(rho_ref_mat.^2, 1));            % 1 x nTrials
                                        for i = 1:S
                                            u_i = reshape(hat_u_all(:,i,:,d_idx), K, nTrials);   % K x nTrials
                                            num_r = sum(rho_ref_mat .* u_i, 1);                  % 1 x nTrials
                                            den_r = sqrt(sum(u_i.^2, 1)) .* rho_ref_norm;        % 1 x nTrials
                                            norm_corr_all(1,i,:,d_idx) = reshape(num_r ./ den_r, 1, 1, nTrials);
                                        end
                                    end

                                    norm_corr_flagged = norm_corr_all < corr_min;
                                    flagged = flagged | norm_corr_flagged;
                                    
                                    %% --- Individual (per-sensor) t0 and alpha estimation ---
                                    mi_sq = reshape(mi_5d.^2, S, 1, nDeployments);
                                    t0_i_hat    = reshape(t(t0_col_idx_per_sensor), 1, S, nTrials, nDeployments);
                                    alpha_i_hat = reshape(c_scalar_at_best ./ (G_scalar_at_best .* mi_sq), 1, S, nTrials, nDeployments);

                                    %% --- Compute trimmed median estimates
                                    k_trim = floor(S/4);
                                    sorted_alpha = sort(alpha_i_hat, 2);
                                    sorted_t0 = sort(t0_i_hat, 2);
                                    alpha_estimates_trimmed = median(sorted_alpha, 2); % median(sorted_alpha(1,k_trim+1:end-k_trim,:,:), 2);
                                    t0_estimates_trimmed = median(sorted_t0, 2); % median(sorted_t0(1,k_trim+1:end-k_trim,:,:), 2);

                                    %% --- MAD Test
                                    med_alpha = median(alpha_i_hat, 2);   % 1 x S x nTrials x nDeployments
                                    mad_alpha = mad(alpha_i_hat, 1, 2);   % median absolute deviation across sensors
                                    med_t0 = median(t0_i_hat, 2);
                                    mad_t0 = mad(t0_i_hat, 1, 2);

                                    tol_pct = 0.15;
                                    alpha_tol = alpha_true * tol_pct;
                                    t0_tol = t0_true * tol_pct;
                                    t0_idx_tol = ceil(t0_true * tol_pct/dt);
                                    
                                    %%% With alpha
                                    alpha_delta_k = 5;
                                    alpha_z = (alpha_i_hat - med_alpha);
                                    alpha_flagged = abs(alpha_z) > max(alpha_tol, alpha_delta_k*mad_alpha);
                                    flagged = flagged | alpha_flagged;

                                    %%% With t0
                                    t0_delta_k = 5;
                                    t0_z = t0_i_hat - med_t0;
                                    t0_flagged = abs(t0_z) > max(t0_tol, t0_delta_k*mad_t0);
                                    flagged = flagged | t0_flagged;

                                    %% Nominal estimation using multiple receive antennas
                                    % The t0 search and the alpha estimate are both linear in the
                                    % per-branch filtered data, so build v = sum_m b_m*z_m*Qn_m and
                                    % Q_sum = sum_m b_m^2*Qn_m once and get every estimate from those.
                                    W_bcast = reshape(W, Mtot, Mtot, 1, nDeployments);
                                    z    = pagemtimes(W_bcast, y_R);      % 2M x K x nTrials x nDeployments
                                    z_bs = pagemtimes(W_bcast, y_R_bs);

                                    mu = m.^2;
                                    WG_Rmu = reshape(pagemtimes(W, pagemtimes(G_R, reshape(mu, S, 1, nDeployments))), Mtot, nDeployments);

                                    v     = zeros(nTrials, K, nDeployments);
                                    v_bs  = zeros(nTrials, K, nDeployments);
                                    Q_sum = zeros(K, K, nDeployments);
                                    for m_idx = 1:Mtot
                                        lambda_m = reshape(lambda_vals(m_idx,:), 1, 1, 1, nDeployments);
                                        [~, Qn_m] = get_time_domain(Hm_arr(omega, lambda_m), dt, nDeployments);   % K x K x D
                                        b_m = reshape(WG_Rmu(m_idx,:), 1, 1, nDeployments);

                                        v     = v     + b_m .* pagemtimes(permute(z(m_idx,:,:,:),    [3 2 4 1]), Qn_m);
                                        v_bs  = v_bs  + b_m .* pagemtimes(permute(z_bs(m_idx,:,:,:), [3 2 4 1]), Qn_m);
                                        Q_sum = Q_sum + b_m.^2 .* Qn_m;
                                    end

                                    mf_t0 = pagemtimes(v, D_template);       mf_t0_bs = pagemtimes(v_bs, D_template);
                                    [~, I]    = max(mf_t0(:, offset_idx:end, :),    [], 2);
                                    [~, I_bs] = max(mf_t0_bs(:, offset_idx:end, :), [], 2);

                                    t0_col    = I    + offset_idx - 1;
                                    t0_estimates_for_plot = reshape(t(t0_col),                 1, 1, nTrials, nDeployments);
                                    t0_estimates_bs       = reshape(t(I_bs + offset_idx - 1),  1, 1, nTrials, nDeployments);

                                    % Alpha given t0. Joint (default): each estimator's own t0_hat, which lies on the
                                    % grid, so rho(t - t0_hat) is a column of D_template. Otherwise: the true t0.
                                    if alpha_uses_true_t0
                                        Y_true = dt*ifft(fft(sensor_signal(t - t0_true, Tp, norm_fact), N, 1) .* mfTemplateFFT_raw);
                                        Rss    = repmat(Y_true(round(Tp/dt)+1 : round(Tp/dt)+K).', nTrials, 1, nDeployments);   % nTrials x K x D
                                        Rss_bs = Rss;
                                    else
                                        Rss    = permute(reshape(D_template(:, t0_col(:)),                K, nTrials, nDeployments), [2 1 3]);
                                        Rss_bs = permute(reshape(D_template(:, I_bs(:) + offset_idx - 1), K, nTrials, nDeployments), [2 1 3]);
                                    end
                                    denom    = dt*dt * sum(pagemtimes(Rss,    Q_sum) .* Rss,    2);   % nTrials x 1 x D
                                    denom_bs = dt*dt * sum(pagemtimes(Rss_bs, Q_sum) .* Rss_bs, 2);
                                    alpha_estimates    = reshape(dt*dt*sum(v    .* Rss,    2) ./ denom,    1, 1, nTrials, nDeployments);
                                    alpha_estimates_bs = reshape(dt*dt*sum(v_bs .* Rss_bs, 2) ./ denom_bs, 1, 1, nTrials, nDeployments);

                                    % % Alpha: t0_hat is on the grid, so Rss(t, t0_hat) is a column of D_template (no FFT).
                                    % Rss = permute(reshape(D_template(:, t0_col(:)), K, nTrials, nDeployments), [2 1 3]);   % nTrials x K x D
                                    % denom = dt*dt * sum(pagemtimes(Rss, Q_sum) .* Rss, 2);    % nTrials x 1 x D
                                    % alpha_estimates    = reshape(dt*dt*sum(v    .* Rss, 2) ./ denom, 1, 1, nTrials, nDeployments);
                                    % alpha_estimates_bs = reshape(dt*dt*sum(v_bs .* Rss, 2) ./ denom, 1, 1, nTrials, nDeployments);

                                    % --- Baseline: full S-sensor array, as if the attacker never
                                    % hijacked anyone
                                    alpha_estimates_baseline = zeros(size(alpha_estimates));
                                    t0_estimates_baseline    = zeros(size(t0_estimates_for_plot));
                                    for d_idx = 1:nDeployments
                                        geom_baseline = build_null_geometry([], mi_5d, g_tilde, gamma_w, gamma_n, Hm_arr, omega, dt, K, S, d_idx);
                                        [alpha_b, t0_b, ~] = estimate_with_geometry(geom_baseline, y_no_attack, K, N, Tp, ...
                                            norm_fact, t, t0_true, mfTemplateFFT_raw, D_template, 1:nTrials, d_idx, dt, offset_idx, alpha_uses_true_t0);
                                        alpha_estimates_baseline(1,1,:,d_idx) = alpha_b;
                                        t0_estimates_baseline(1,1,:,d_idx) = t0_b;
                                    end

                                    % Per-trial flagged sensor SET (face value -- no validation).
                                    % May be empty, a single sensor, or multiple sensors if more
                                    % than one clears threshold in a given trial.
                                    excluded_sets = cell(nTrials, nDeployments);
                                    for d_idx = 1:nDeployments
                                        for tr = 1:nTrials
                                            flags_this = squeeze(flagged(1,:,tr,d_idx));   % 1 x S logical
                                            excluded_sets{tr,d_idx} = find(flags_this);     % row vec, possibly empty
                                        end
                                    end

                                    % Group trials within each deployment by their EXACT flagged
                                    % set (usually just a couple of distinct groups given how
                                    % separated T_i is), and null/re-estimate once per group --
                                    % avoids one function call per individual trial.
                                    if attacker_enabled
                                        alpha_estimates_undefended = alpha_estimates;
                                        t0_estimates_undefended    = t0_estimates_for_plot;

                                        % --- Oracle ceiling: null the TRUE attacker_idx directly
                                        alpha_estimates_oracle = zeros(size(alpha_estimates));
                                        t0_estimates_oracle    = zeros(size(t0_estimates_for_plot));
                                        for d_idx = 1:nDeployments
                                            geom_oracle = build_null_geometry(attacker_idx, mi_5d, g_tilde, gamma_w, gamma_n, Hm_arr, omega, dt, K, S, d_idx);
                                            [alpha_o, t0_o, ~] = estimate_with_geometry(geom_oracle, y, K, N, Tp, ...
                                                norm_fact, t, t0_true, mfTemplateFFT_raw, D_template, 1:nTrials, d_idx, dt, offset_idx, alpha_uses_true_t0);
                                            alpha_estimates_oracle(1,1,:,d_idx) = alpha_o;
                                            t0_estimates_oracle(1,1,:,d_idx) = t0_o;
                                        end
                                    end

                                    alpha_final = alpha_estimates;
                                    t0_final = t0_estimates_for_plot;

                                    % --- Face value validation
                                    num_geometries_built = 0;
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
                                            tp_total = tp_total + numel(intersect(A, attacker_idx)) * n_members;
                                            fp_total = fp_total + numel(setdiff(A, attacker_idx))   * n_members;
                                            fn_total = fn_total + numel(setdiff(attacker_idx, A))   * n_members;

                                            if isempty(A)
                                                continue   % nothing flagged -- keep undefended estimate
                                            end

                                            geom = build_null_geometry(A, mi_5d, g_tilde, gamma_w, gamma_n, Hm_arr, omega, dt, K, S, d_idx);
                                            num_geometries_built = num_geometries_built + 1;

                                            [alpha_A, t0_A, ~] = estimate_with_geometry(geom, y, K, N, Tp, ...
                                                norm_fact, t, t0_true, mfTemplateFFT_raw, D_template, trial_members, d_idx, dt, offset_idx, alpha_uses_true_t0);

                                            alpha_final(1,1,trial_members,d_idx) = alpha_A;
                                            t0_final(1,1,trial_members,d_idx) = t0_A;
                                        end
                                    end

                                    log_msg(verbosity_level, 4, 'Face-value nulling: %d distinct null-geometries built across all deployments', num_geometries_built);
                                    log_msg(verbosity_level, 4, 'Detection accuracy: TP=%d, FP=%d, FN=%d', tp_total, fp_total, fn_total);

                                    honest_idx_diag = setdiff(1:S, attacker_idx);
                                    
                                    % --- FP rate lookup: one struct/map holding every test's honest-sensor FP
                                    % rate, keyed by name. Add a test with one line; call fp_stats('name') or
                                    % dump the whole thing with fp_report(fp_stats) any time you want a look.
                                    fp_stats = containers.Map('KeyType','char','ValueType','double');
                                    fp_stats('alpha_i')          = mean(alpha_flagged(1,honest_idx_diag,:,:), 'all');
                                    fp_stats('t0_i')             = mean(t0_flagged(1,honest_idx_diag,:,:), 'all');
                                    fp_stats('normcorr')         = mean(norm_corr_flagged(1,honest_idx_diag,:,:), 'all');
                                    fp_stats('combined')         = mean(flagged(1,honest_idx_diag,:,:), 'all');
                                    
                                    fp_report(fp_stats, verbosity_level);

                                    alpha_estimates = alpha_final;
                                    t0_estimates_for_plot = t0_final;
                                end
                                
                                slot = {strat_idx, agent_db_idx, channel_db_idx, ':', scheme_idx, pivot_idx, save_dim, ':'};
                                rho_empirical_var(slot{:}) = [reshape(var(alpha_estimates,0,3),1,[]); reshape(var(t0_estimates_for_plot,0,3),1,[])];
                                rho_empirical_mse(slot{:}) = alpha_t0_mse(alpha_estimates, t0_estimates_for_plot, alpha_true, t0_true);

                                if attacker_enabled && use_W
                                    rho_empirical_mse_baseline(slot{:})   = alpha_t0_mse(alpha_estimates_baseline,   t0_estimates_baseline,   alpha_true, t0_true);
                                    rho_empirical_mse_undefended(slot{:}) = alpha_t0_mse(alpha_estimates_undefended, t0_estimates_undefended, alpha_true, t0_true);
                                    rho_empirical_mse_oracle(slot{:})     = alpha_t0_mse(alpha_estimates_oracle,     t0_estimates_oracle,     alpha_true, t0_true);
                                    rho_empirical_mse_bs(slot{:})         = alpha_t0_mse(alpha_estimates_bs,         t0_estimates_bs,         alpha_true, t0_true);
                                    rho_empirical_mse_trimmed(slot{:})    = alpha_t0_mse(alpha_estimates_trimmed,    t0_estimates_trimmed,    alpha_true, t0_true);
                                end

                                if use_W == false
                                    rho_crlb(strat_idx,agent_db_idx,channel_db_idx,1,scheme_idx,pivot_idx,save_dim,:) = alpha_objective_array(1, mi_4d, ci);
                                    rho_crlb(strat_idx,agent_db_idx,channel_db_idx,2,scheme_idx,pivot_idx,save_dim,:) = t0_objective_array(1, mi_4d, ci, alpha_true);
                                else
                                    b_all   = WG_Rmu;         % Mtot x D
                                    lam_all = lambda_vals;    % Mtot x D
                                    alpha_term = integral(@(w) sum(b_all.^2 .*          mag_sqr_S0(w).^2 ./ (lam_all.*mag_sqr_S0(w)+1), 1), -Inf, Inf, 'ArrayValued', true);
                                    t0_term    = integral(@(w) sum(b_all.^2 .* (w.^2) .* mag_sqr_S0(w).^2 ./ (lam_all.*mag_sqr_S0(w)+1), 1), -Inf, Inf, 'ArrayValued', true);

                                    rho_crlb(strat_idx,agent_db_idx,channel_db_idx,1,scheme_idx,pivot_idx,save_dim,:) = 2*pi ./ alpha_term;
                                    rho_crlb(strat_idx,agent_db_idx,channel_db_idx,2,scheme_idx,pivot_idx,save_dim,:) = 2*pi ./ (alpha_true^2 .* t0_term);

                                    if attacker_enabled
                                        alpha_term_oracle = integral(@(w) sum(WG_Rmu_oracle.^2 .*          mag_sqr_S0(w).^2 ./ (lambda_vals_oracle.*mag_sqr_S0(w)+1), 1), -Inf, Inf, 'ArrayValued', true);
                                        t0_term_oracle    = integral(@(w) sum(WG_Rmu_oracle.^2 .* (w.^2) .* mag_sqr_S0(w).^2 ./ (lambda_vals_oracle.*mag_sqr_S0(w)+1), 1), -Inf, Inf, 'ArrayValued', true);

                                        rho_crlb_oracle(strat_idx,agent_db_idx,channel_db_idx,1,scheme_idx,pivot_idx,save_dim,:) = 2*pi ./ alpha_term_oracle;
                                        rho_crlb_oracle(strat_idx,agent_db_idx,channel_db_idx,2,scheme_idx,pivot_idx,save_dim,:) = 2*pi ./ (alpha_true^2 .* t0_term_oracle);
                                    end
                                end
                            end
                            total_est_time = toc(total_est_start);
                            log_msg(verbosity_level, 3, 'Total ML estimation time: %dm %.1fs', floor(total_est_time/60), mod(total_est_time,60));
                        
                            total_channel_db_time = toc(channel_db_start);
                            log_msg(verbosity_level, 3, 'Runtime for channel SNR = %g dB: %dm %.1fs\n', ...
                                channel_db_values(scheme_idx,channel_db_idx,agent_db_idx), floor(total_channel_db_time/60), mod(total_channel_db_time,60));
                        end
                    end
                end
            end
        end
    end
    
    % Stop run timer.
    runtime = toc(loopTic);
    % Display total run time.
    log_msg(verbosity_level, 1, '%s took %dm %.1fs', experiment, floor(runtime/60), mod(runtime,60));
    
    %% Compute averages across deployments.
    if nDeployments == 1
        avg_dep_var = rho_empirical_var;
        avg_dep_mse = rho_empirical_mse;
        avg_dep_crlb = rho_crlb;
        avg_dep_bias = rho_empirical_bias;
        avg_dep_mse_undefended = rho_empirical_mse_undefended;
        avg_dep_mse_oracle = rho_empirical_mse_oracle;   % (or mean(...) branch, matching the existing if/else)
        avg_dep_mse_baseline = rho_empirical_mse_baseline;
        avg_dep_mse_trimmed = rho_empirical_mse_trimmed;
        avg_dep_mse_bs = rho_empirical_mse_bs;
        avg_dep_crlb_oracle = rho_crlb_oracle;
    else
        avg_dep_var = mean(rho_empirical_var,length(size(rho_empirical_var)));
        avg_dep_mse = mean(rho_empirical_mse,length(size(rho_empirical_mse)));
        avg_dep_crlb = mean(rho_crlb,length(size(rho_crlb)));
        avg_dep_bias = mean(rho_empirical_bias,length(size(rho_empirical_bias)));
        avg_dep_mse_undefended = mean(rho_empirical_mse_undefended,length(size(rho_empirical_mse_undefended)));
        avg_dep_mse_oracle = mean(rho_empirical_mse_oracle,length(size(rho_empirical_mse_oracle)));   % (or mean(...) branch, matching the existing if/else)
        avg_dep_mse_baseline = mean(rho_empirical_mse_baseline,length(size(rho_empirical_mse_baseline)));
        avg_dep_mse_trimmed = mean(rho_empirical_mse_trimmed,length(size(rho_empirical_mse_trimmed)));
        avg_dep_mse_bs = mean(rho_empirical_mse_bs,length(size(rho_empirical_mse_bs)));
        avg_dep_crlb_oracle = mean(rho_crlb_oracle,length(size(rho_crlb_oracle)));
    end

    txt_list = ["undefended", "oracle", "baseline", "trimmed", "bs"];
    save_images = false;

    files_to_save = ["avg_dep_mse_" + txt_list, "avg_dep_var", "avg_dep_mse", "avg_dep_bias", "avg_dep_crlb", "avg_dep_crlb_oracle",...
        "agent_db_values", "channel_snr", "selected_schemes", "dropout_vals", "sensor_vals", "channel_db_values",...
        "experiment", "sensor_dimension", "save_images", "K", "nDeployments", "constraints", "transforms", "coeff_types", "approaches",...
        "iter_arr", "iter_length", "all_strats", "num_strats", "save_images", "alpha_uses_true_t0",...
        "T0", "Tp", "alpha_true", "t0_true", "pivot_vals", "Mtot", "N", "S", "S_max","attacks", "attacker_enabled"];
    
    save(experiment + "_results.mat", files_to_save{:})
    
    %% Plotting placeholder
    plot_attacks
end