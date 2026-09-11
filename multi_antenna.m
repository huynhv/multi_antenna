% 9/3/2026 LEFT OFF HERE: atp it seems like jamming and flipping are
% solved, the more open-ended attacks are the arbitrary time shift and
% amplitude scaling

% Clear all variables and close all existing figures.
clearvars
close all

addpath('.\functions')
verbosity_level = 5;   % 1 = major sections only, 2 = + sub-steps/channel-SNR, 3 = + diagnostics

% Set seed for reproducibility
seed = 2048;
rng(seed);

nDeployments = 10;
% Set the number of measurements per sensor, starting from 0
K = 301;
% Set number of trials.
nTrials = 500;

% Set number of sensors
sensor_vals = 5; % [5,10]
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
t0_true = 1.5;

% Set the number of antennas
num_antennas = 10;
Mtot = 2*num_antennas;

%%% Always generate w, n, and channel gain for maximum S, then use subsets for different S trials
S_max = max(sensor_vals);
% Generate sensor noise with dimensions: K x S x trials x deployments.
all_w = randn(K,S_max,1,nTrials,nDeployments);
% Generate server noise.
n = (randn(K,1,num_antennas,nTrials,nDeployments) + 1i*randn(K,1,num_antennas,nTrials,nDeployments));
% Generate channel gains.
all_gi = (randn(1,S_max,num_antennas,1,nDeployments) + 1i*randn(1,S_max,num_antennas,1,nDeployments))/sqrt(2);

% Define range for distribution of mi and ti (tau_i)
% max_mi = 1;
% min_mi = 0.5;

max_ti = 1;
min_ti = 0.5; 

%% Generate mi and ti
% Distance-first: ti (= di, v=1) sampled directly since max_ti is the true
% structural constraint (bounded by the fixed observation window). mi is
% DERIVED from di via a physical path-loss + shadowing model, so it remains
% an independent, non-circular quantity for the geometric defense to test
% against later (g(d_i) vs. the sensor's implied amplitude).
path_loss_exp = 1;
shadow_sigma_dB = 1;

m_ref = 1;              % reference amplitude at distance min_ti
c_gain = m_ref * min_ti^path_loss_exp;   % anchors: unshadowed mi = m_ref at d = min_ti

all_di = sort(unifrnd(min_ti, max_ti, 1, S_max, 1, 1, nDeployments), 2, 'ascend');
all_ti = all_di;   % v = 1, per earlier discussion

shadowing = 10.^(shadow_sigma_dB/20 * randn(1,S_max,1,1,nDeployments));
all_mi = c_gain ./ (all_di.^path_loss_exp) .* shadowing;

% Diagnostic: confirm mi's mean stays representative (mean/median should be
% reasonably close; large divergence signals a heavy tail that would make
% E_mi_sqr, and hence gamma_n/Ki, poorly calibrated -- same issue found
% earlier when mi was fully unbounded).
log_msg(verbosity_level, 3, 'mi: mean=%.4f, median=%.4f, std=%.4f', ...
    mean(all_mi(:)), median(all_mi(:)), std(all_mi(:)));
log_msg(verbosity_level, 3, 'E[mi^2] via mean vs. median^2: %.4f vs %.4f', ...
    mean(all_mi(:).^2), median(all_mi(:))^2);

% Deterministic path-loss spread vs. shadowing -- shadowing should be a
% modest perturbation, not large enough to make near/far reordering routine.
det_spread_dB = 20*log10( (c_gain/min_ti^path_loss_exp) / (c_gain/max_ti^path_loss_exp) );
log_msg(verbosity_level, 3, 'Deterministic path-loss spread: %.1f dB | Shadowing std dev: %.1f dB', ...
    det_spread_dB, shadow_sigma_dB);

%%
selected_schemes = "EPC";

% Define Rayleigh distribution parameters.
rayleigh_factor = 1/sqrt(2);
E_mag_g_sqr = 2*rayleigh_factor^2;

% Define expected value of mi^2.
% E_mi_sqr = (1/12) * (max_mi-min_mi)^2 + 0.5*(min_mi+max_mi);

E_mi_sqr = mean(all_mi(:).^2);

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

%% Byzantine attacker configuration
% Attacker(s) hijack existing sensor slots and transmit unstructured (random)
% noise in place of the honest matched-filter output. Adaptable to multiple
% simultaneous attackers by extending attacker_idx.

% Generates all_attacker_noise per attack_type -- same shape convention the
% rest of the pipeline already expects (K x S_max x 1 x nTrials x nDeployments),
% so gamma_w_attacker/scaled_attacker_noise/the injection loop below need no
% changes. Deterministic attacks (dc/square/sawtooth/sinusoid) draw ONE
% independent set of parameters per potential attacker slot, reused across
% every trial; only "noise" is redrawn per trial.
attacker_enabled = true;
attacker_idx     = [1];      % sensor index/indices (within 1:S) that are compromised
attacker_db      = -10;      % attacker "SNR", same convention as agent_db (see below)
attacker_seed = 42;

use_greedy_validation = false;   % true: Lambda-validated greedy; false: face-value nulling
peak_same_location = true;
attack_type = "noise";   % "noise" | "dc" | "square" | "sawtooth" | "sinusoid" | "peak" | "flip" | "replay" | "amplitude"
attacker_beta = 10;
rng(attacker_seed);

f_nominal = Bee/(2*pi);   % = 1/(2*Tp), reference frequency near the honest passband

switch attack_type
    case "noise"
        all_attacker_noise = randn(K,S_max,1,nTrials,nDeployments);

    case "dc"
        dc_sign = sign(randn(1,S_max));                    % 1 x S_max, one sign per slot
        waveform = repmat(dc_sign, K, 1);                  % K x S_max
        all_attacker_noise = repmat(reshape(waveform,K,S_max,1,1,1), 1,1,1,nTrials,nDeployments);

    case "square"
        f_c   = f_nominal * (0.5 + rand(1,S_max));         % 1 x S_max
        phase = 2*pi*rand(1,S_max);
        waveform = sign(sin(2*pi*t.*f_c + phase));         % K x S_max
        all_attacker_noise = repmat(reshape(waveform,K,S_max,1,1,1), 1,1,1,nTrials,nDeployments);

    case "sawtooth"
        f_c        = f_nominal * (0.5 + rand(1,S_max));
        phase_frac = rand(1,S_max);
        phi = t.*f_c + phase_frac;
        waveform = 2*(phi - floor(phi + 0.5));             % K x S_max, in [-1,1)
        all_attacker_noise = repmat(reshape(waveform,K,S_max,1,1,1), 1,1,1,nTrials,nDeployments);

    case "sinusoid"
        f_c   = f_nominal * (0.5 + rand(1,S_max));
        phase = 2*pi*rand(1,S_max);
        waveform = sin(2*pi*t.*f_c + phase);               % K x S_max
        all_attacker_noise = repmat(reshape(waveform,K,S_max,1,1,1), 1,1,1,nTrials,nDeployments);

    case "peak"
        % Concentrated single-sample spike. Location toggle:
        %   peak_same_location = true  -> all attacker slots use ONE shared
        %       draw (matches the manuscript's single-attacker worst case)
        %   peak_same_location = false -> each slot independently randomized
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

    case "replay"
        % In-subspace / coherent replay: beta * rho(t - tau), tau drawn from
        % the feasible delay window T = [0, t0_max] -- the SAME window the
        % honest subspace itself is built from. By construction this is
        % indistinguishable from an honest signal to EITHER defense (basis
        % suppression or T_i) -- included as the known negative control /
        % worst case, not something either defense is expected to catch.
        t0_grid_full = repmat(reshape(t,1,1,[]), 1,1,1,nDeployments);
        [~, R00_full] = mf_integral_fft(sensor_signal(t-t0_grid_full,Tp,norm_fact), sensor_signal(t,Tp,norm_fact), 1, 1, K, dt, Tp);
        rho_dictionary = reshape(R00_full(:,:,:,1), K, K);   % rho(t - t_j), same construction as D_template

        feasible_grid = t(t <= t0_max);
        tau_per_slot = feasible_grid(randi(numel(feasible_grid), 1, S_max));
        % tau_per_slot = 0.1 * ones(1, S_max); 

        waveform = zeros(K, S_max);
        for a = 1:S_max
            [~, col_idx] = min(abs(t - tau_per_slot(a)));
            waveform(:,a) = rho_dictionary(:, col_idx);
        end
        all_attacker_noise = repmat(reshape(waveform,K,S_max,1,1,1), 1,1,1,nTrials,nDeployments);

    case "flip"
        % Handled entirely at injection time below (needs the actual per-trial
        % honest u_i(t), which doesn't exist yet at this point in the script).
        % No independent waveform to generate or energy-normalize here.
        all_attacker_noise = zeros(K,S_max,1,nTrials,nDeployments);   % unused placeholder

    case "amplitude"
        % Same reasoning as "flip" -- this attack SCALES the actual honest
        % u_i(t) rather than transmitting an independent waveform, so it must
        % be handled at injection time below (needs the real per-trial honest
        % signal, which doesn't exist yet here). No independent waveform to
        % generate or energy-normalize.
        all_attacker_noise = zeros(K,S_max,1,nTrials,nDeployments);   % unused placeholder

    otherwise
        error("Unknown attack_type: %s", attack_type);
end

%% Normalize EVERY attack type to unit total (dt-weighted) energy per slot,
% BEFORE the shared attacker_db/gamma_w_attacker scaling is applied below.
% Without this, concentrated waveforms (e.g. "peak") end up with far less
% TOTAL energy than spread-out ones (noise/square/etc.) under the same
% per-sample scale factor -- this keeps the energy convention genuinely
% consistent across every attack type, matching what all attacks are
% supposed to share.
if attack_type ~= "flip"
    waveform_energy = squeeze(sum(all_attacker_noise(:,:,1,1,1).^2, 1)) * dt;   % 1 x S_max
    norm_factor = reshape(1./sqrt(waveform_energy), 1, S_max, 1, 1, 1);
    all_attacker_noise = all_attacker_noise .* norm_factor;
end

rng(seed);

w_psd_constant = 1;
n_psd_constant = 1;
N = 2*K - 1;

%% Subspace-projection residual test: honest-signal basis Phi and its
% orthogonal complement Psi. Built ONCE, offline -- depends only on the known
% pulse shape and the K-point time grid, never on S, gains, or trial data.
t0_grid_top = repmat(reshape(t,1,1,[]), 1,1,1,nDeployments);
[~, R00_full] = mf_integral_fft(sensor_signal(t-t0_grid_top,Tp,norm_fact), sensor_signal(t,Tp,norm_fact), 1, 1, K, dt, Tp);
D_template = reshape(R00_full(:,:,:,1), K, K);   % K x K, [D_template]_{k,l} = rho(t_k - t_l)

[U_D, Sigma_D, ~] = svd(D_template, 'econ');
sv_energy = cumsum(diag(Sigma_D).^2) / sum(diag(Sigma_D).^2);
d_sub = find(sv_energy >= 0.999, 1, 'first');
Phi = U_D(:, 1:d_sub);           % K x d_sub: honest-signal subspace
Psi = U_D(:, d_sub+1:end);       % K x (K-d_sub): its orthogonal complement -- free byproduct of the same SVD

% Discretized autocorrelation restricted to Psi-coordinates -- the known
% "shape" of colored sensor noise projected into the residual space, needed
% below to build each sensor's noise covariance K_i.
Rss_Psi = Psi.' * D_template * Psi;   % (K-d_sub) x (K-d_sub)

% Eigendecompose Rss_Psi ONCE, offline -- K_i = a*Rss_Psi + b*I always shares
% these eigenvectors for ANY a,b, so K_i^-1 reduces to a cheap per-eigenvalue
% reweight instead of a fresh (K-d_sub)^3 matrix inverse every time K_i is needed.
[V_Rss, Lam_Rss] = eig((Rss_Psi+Rss_Psi.')/2);
lam_Rss = diag(Lam_Rss);   % (K-d_sub) x 1

% Eigenvalues of D_template itself (already computed via U_D/Sigma_D above) --
% needed below to whiten against a single candidate direction rho(.-t0_hat)
% rather than the full Phi/Psi split, for the coherent-replay (T_i^align) test.
lam_D = diag(Sigma_D);   % K x 1

log_msg(verbosity_level, 1, 'Subspace-projection residual test: d = %d, K-d = %d', d_sub, K-d_sub);

%% Define Anonymous Expressions
mag_sqr_S0_internal = @(w) (norm_fact^2) * (2*Bee.^2 .* (1 + cos(w.*pi./Bee)) ) ./ (w.^2 - Bee.^2).^2;
mag_sqr_S0 = @(w) zeroIfnan(mag_sqr_S0_internal(w)) + (w == Bee | w == -Bee) .* (mag_sqr_S0_internal(w-epsilon) + mag_sqr_S0_internal(w+epsilon))./2;

%%% Expressions for original manuscript
V_arr_internal = @(w, ai, mi, ci) scaled_w_psd_constant .* sum(ci.^2 .* ai.^2 .* mi.^2, 2) .* mag_sqr_S0_internal(w) + scaled_n_psd_constant/2;
V_arr = @(w, ai, mi, ci) zeroIfnan(V_arr_internal(w, ai, mi, ci)) + (w == Bee | w == -Bee) .* (V_arr_internal(w-epsilon, ai, mi, ci) + V_arr_internal(w+epsilon, ai, mi, ci))./2;

H_arr_internal = @(w, ai, mi, ci) (1 + epsilon) ./ (V_arr_internal(w, ai, mi, ci) + epsilon);
H_arr = @(w, ai, mi, ci) zeroIfnan(H_arr_internal(w, ai, mi, ci)) + (w == Bee | w == -Bee) .* (H_arr_internal(w-epsilon, ai, mi, ci) + H_arr_internal(w+epsilon, ai, mi, ci))./2;

alpha_objective_array = @(ai, mi, ci) (2*pi) ./( ((sum(ci .* ai .* mi.^2)).^2) .* integral(@(w) (mag_sqr_S0(w)).^2 ./ V_arr(w, ai, mi, ci), -Inf, Inf, 'ArrayValued', true));
t0_objective_array = @(ai, mi, ci, arb_alpha) (2*pi) ./( (arb_alpha^2 .* (sum(ci .* ai .* mi.^2)).^2) .* integral(@(w) (w .* mag_sqr_S0(w)).^2 ./ V_arr(w,ai, mi, ci), -Inf, Inf, 'ArrayValued', true));

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
    channel_snr = linspace(0,15,5);
    rho_vals = 0.5;

    pivot_vals = rho_vals;
    
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
    rho_crlb = zeros(num_strats,length(agent_db_values),size(channel_db_values,2),2,length(selected_schemes),length(rho_vals),sensor_dimension,nDeployments);
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
        
        for pivot_idx = 1:length(pivot_vals)
        
            rho = rho_vals(pivot_idx);
            obj_func = "crlb sum";

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

                    % Capture the PRISTINE, never-attacked ui -- used below to build the
                    % "baseline" (full S-sensor array, attacker never hijacked anyone)
                    % comparison. NOT related to the "single-antenna baseline" terminology
                    % used elsewhere in this script -- kept as ui_no_attack to avoid confusion.
                    if attacker_enabled
                        ui_no_attack = ui;
                    end

                    %% --- Byzantine attacker: overwrite compromised sensor(s)' transmitted
                    % signal with unstructured noise, in place of their honest ui(t). ---
                    if attacker_enabled
                        active_attackers = attacker_idx(attacker_idx <= S);
                        if attack_type == "flip"
                            % Sign-flip: -alpha * m_i^2 * rho(t-t0) - w_i(t). Still a
                            % LINEAR member of the honest subspace S (span includes
                            % negative scalings) -- provably undefeatable by either
                            % defense, same category as "replay". Included as a
                            % negative control, not a genuine test of detection.
                            for a = active_attackers
                                ui(:, a, :, :, :) = -ui(:, a, :, :, :);
                            end
                        elseif attack_type == "amplitude"
                            % General amplitude lie: beta * alpha * m_i^2 * rho(t-t0) +
                            % beta * w_i(t) -- CORRECT t0, CORRECT shape, WRONG magnitude.
                            % Passes T_i (correct shape), sign-flip (positive if beta>0),
                            % and T_i^align (correct timing) -- undefeatable by any
                            % content-only test. Target case for the geometric/position
                            % defense, not yet implemented.
                            for a = active_attackers
                                ui(:, a, :, :, :) = attacker_beta * ui(:, a, :, :, :);
                            end
                        else
                            for a = active_attackers
                                ui(:, a, :, :, :) = scaled_attacker_noise(:, a, :, :, :);
                            end
                        end
                    end
                    % [~, ui_noise] = mf_integral_fft(scaled_w(:,1:S,:,:,:), mi_5d .* sensor_signal(t, Tp, norm_fact), 1, 1, K, dt, Tp);
    
                    use_W = true;

                    % Define search space for t0
                    t0_search_space = repmat(reshape(t,1,1,[]), 1,1,1,nDeployments);
                    if ~use_W
                        [~, ui_t0_search] = mf_integral_fft(mi_4d .* sensor_signal(t-t0_search_space, Tp, norm_fact), mi_4d .* sensor_signal(t, Tp, norm_fact), 1, 1, K, dt, Tp);
                        mfTemplate = mi_4d .* sensor_signal(t, Tp, norm_fact);
                        mfTemplateFFT = fft(mfTemplate, N, 1);
                    else
                        [~, R00_t0_search] = mf_integral_fft(sensor_signal(t-t0_search_space,Tp,norm_fact), sensor_signal(t,Tp,norm_fact), 1, 1, K, dt, Tp);  % K x 1 x K x D
                        mfTemplateFFT_raw = fft(sensor_signal(t, Tp, norm_fact), N, 1);
                    end

                    if scheme == "EPC"
                        % Apply exact phase compensation
                        ref_ant_idx = 1;
                        g_ref = gi(:, :, ref_ant_idx, :, :);
                        phase_comp = exp(-1j * angle(g_ref));
                        g_tilde = gi .* phase_comp;
                    end
                        
                    out_cws = sum(g_tilde .* ui, 2);
                    % out_cws_noise = sum(g_tilde .* ui_noise, 2);
                    if attacker_enabled
                        out_cws_no_attack = sum(g_tilde .* ui_no_attack, 2);
                    end

                    %% Compute values for multi-antenna
                    m = reshape(mi_5d, S, nDeployments);           % S x nDeployments
                    g = reshape(g_tilde, S, num_antennas, nDeployments);       % S x num_antennas x nDeployments (complex, per-antenna compensated gain)

                    G_R = cat(2, real(g), imag(g));
                    G_R = permute(G_R, [2,1,3]);   % Mtot x S x nDeployments (channel-independent)

                    combined_noise_psd = reshape(m.^2 .* gamma_w, S, 1, nDeployments);
                    Dmat = eye(S) .* combined_noise_psd;   % S x S x nDeployments

                    % ---- channel-INDEPENDENT: spatial kernel + eigenvectors (once per deployment) ----

                    % --- Hermitian part: A = E{v v^H} spatial kernel (M x M x nDeployments) ---
                    A = pagemtimes(pagemtimes(pagetranspose(g), Dmat), conj(g));
                    A = (A + pagectranspose(A)) / 2;

                    % --- Pseudo-covariance part: Atilde = E{v v^T} (no conjugate) ---
                    Atilde = pagemtimes(pagemtimes(pagetranspose(g), Dmat), g);
                    Atilde = (Atilde + pagetranspose(Atilde)) / 2;

                    B11 =  0.5 * real(A + Atilde);
                    B12 = -0.5 * imag(A - Atilde);
                    B21 =  0.5 * imag(A + Atilde);
                    B22 =  0.5 * real(A - Atilde);

                    B = cat(1, cat(2, B11, B12), cat(2, B21, B22));   % 2M x 2M x nDeployments
                    B = (B + pagetranspose(B)) / 2;

                    % Eigen-decompose B ITSELF (not Bprime).  Because Sigma_n_R is a scalar multiple
                    % of I, Bprime = (2/gamma_n)*B, so the eigenVECTORS are those of B and only the
                    % eigenVALUES scale with channel SNR.
                    [U_B, Lam_B] = pageeig(B);
                    for d_idx = 1:nDeployments
                        [lam_sorted, idx] = sort(diag(Lam_B(:,:,d_idx)), 'descend');
                        U_B(:,:,d_idx)      = U_B(:,idx,d_idx);
                        Lam_B(:,:,d_idx)    = diag(lam_sorted);
                    end

                    UB_transpose  = pagetranspose(U_B);               % 2M x 2M x nDeployments  <-- cached "W base"
                    lambda_B_base = zeros(Mtot, nDeployments);        % 2M x nDeployments
                    for d_idx = 1:nDeployments
                        lambda_B_base(:,d_idx) = diag(Lam_B(:,:,d_idx));
                    end

                    % --- Oracle CRLB: rebuild the SAME channel-independent pipeline above,
                    % using ONLY the honest S-1 sensors -- as if the attacker's slot never
                    % existed. No nulling geometry involved; this is the ordinary CRLB
                    % construction, just restricted to a smaller sensor set from the start.
                    if attacker_enabled
                        honest_only = setdiff(1:S, attacker_idx);
                        S_honest = numel(honest_only);

                        m_oracle = m(honest_only,:);                    % S_honest x nDeployments
                        g_oracle = g(honest_only,:,:);                  % S_honest x num_antennas x nDeployments

                        G_R_oracle = cat(2, real(g_oracle), imag(g_oracle));
                        G_R_oracle = permute(G_R_oracle, [2,1,3]);      % Mtot x S_honest x nDeployments

                        combined_noise_psd_oracle = reshape(m_oracle.^2 .* gamma_w, S_honest, 1, nDeployments);
                        Dmat_oracle = eye(S_honest) .* combined_noise_psd_oracle;   % S_honest x S_honest x nDeployments

                        A_oracle = pagemtimes(pagemtimes(pagetranspose(g_oracle), Dmat_oracle), conj(g_oracle));
                        A_oracle = (A_oracle + pagectranspose(A_oracle)) / 2;

                        Atilde_oracle = pagemtimes(pagemtimes(pagetranspose(g_oracle), Dmat_oracle), g_oracle);
                        Atilde_oracle = (Atilde_oracle + pagetranspose(Atilde_oracle)) / 2;

                        B11_o =  0.5 * real(A_oracle + Atilde_oracle);
                        B12_o = -0.5 * imag(A_oracle - Atilde_oracle);
                        B21_o =  0.5 * imag(A_oracle + Atilde_oracle);
                        B22_o =  0.5 * real(A_oracle - Atilde_oracle);

                        B_oracle = cat(1, cat(2, B11_o, B12_o), cat(2, B21_o, B22_o));   % 2M x 2M x nDeployments
                        B_oracle = (B_oracle + pagetranspose(B_oracle)) / 2;

                        [U_B_oracle, Lam_B_oracle] = pageeig(B_oracle);
                        for d_idx = 1:nDeployments
                            [lam_sorted_o, idx_o] = sort(diag(Lam_B_oracle(:,:,d_idx)), 'descend');
                            U_B_oracle(:,:,d_idx)   = U_B_oracle(:,idx_o,d_idx);
                            Lam_B_oracle(:,:,d_idx) = diag(lam_sorted_o);
                        end

                        UB_transpose_oracle  = pagetranspose(U_B_oracle);         % 2M x 2M x nDeployments
                        lambda_B_base_oracle = zeros(Mtot, nDeployments);
                        for d_idx = 1:nDeployments
                            lambda_B_base_oracle(:,d_idx) = diag(Lam_B_oracle(:,:,d_idx));
                        end

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
                            log_msg(verbosity_level, 2, 'Scheme %s | Agent SNR = %d dB | Channel SNR = %g dB | S = %d | Dropout = %d', ...
                                scheme, agent_db_values(agent_db_idx), channel_db_values(scheme_idx,channel_db_idx,agent_db_idx), S, dropout_vals(dropout_idx));

                            if scheme == "EPC"
                                gamma_n = (dt/n_psd_constant) * Ps * E_mi_sqr * E_mag_g_sqr / db2magTen(channel_db_values(scheme_idx,channel_db_idx,agent_db_idx));
                            end

                            %% Set server noise power. Note the factor of 1/dt which is
                            % equivalent to applying an anti-aliasing filter.
                            scaled_n = sqrt( (gamma_n*n_psd_constant/dt) / 2 ) * n;
                            scaled_n_psd_constant = gamma_n * n_psd_constant; % == N0/2

                            % --- Per-sensor known noise covariance for the residual test:
                            % K_i = m_i^2*(Nw/2)*Rss_Psi + (N0/2)/||breve_g_i||^2 * I.
                            % Computed once per channel-SNR point (deterministic given known
                            % parameters -- does NOT depend on trial noise realizations, so
                            % it's reused across all nTrials below). ---
                            Nw_over_2 = scaled_w_psd_constant;   % == Nw/2
                            N0_over_2 = scaled_n_psd_constant / (2*dt);   % == N0/2

                            % K_i = a_i*Rss_Psi + b_i*I shares V_Rss's eigenvectors for ANY
                            % a_i,b_i -- store only the per-sensor, per-deployment eigenvalue
                            % reweighting, not full inverse matrices.
                            
                            % Ki_diag_all = cell(S,1);
                            % for i = 1:S
                            %     a_i = reshape(mi_5d(1,i,1,1,:).^2 * Nw_over_2, 1, nDeployments);
                            %     norm_sq_i = zeros(1,nDeployments);
                            %     for d_idx = 1:nDeployments
                            %         norm_sq_i(d_idx) = sum(breve_g_all{i}(:,d_idx).^2);
                            %     end
                            %     b_i = N0_over_2 ./ norm_sq_i;
                            % 
                            %     Ki_diag_all{i} = 1 ./ (lam_Rss .* a_i + b_i);   % (K-d_sub) x nDeployments
                            % end

                            % Ki_full_eig{i}: SAME a_i,b_i as above, but against D_template's
                            % OWN eigenvalues (lam_D) instead of Rss_Psi's -- gives the full
                            % K-dimensional noise model needed by T_i^align (coherent-replay
                            % test), which whitens against a single direction rho(.-t0_hat),
                            % not the reduced Psi-space T_i already uses.
                            Ki_diag_all = cell(S,1);
                            Ki_full_eig = cell(S,1);
                            for i = 1:S
                                a_i = reshape(mi_5d(1,i,1,1,:).^2 * Nw_over_2, 1, nDeployments);
                                norm_sq_i = zeros(1,nDeployments);
                                for d_idx = 1:nDeployments
                                    norm_sq_i(d_idx) = sum(breve_g_all{i}(:,d_idx).^2);
                                end
                                b_i = N0_over_2 ./ norm_sq_i;

                                Ki_diag_all{i} = 1 ./ (lam_Rss .* a_i + b_i);   % (K-d_sub) x nDeployments
                                Ki_full_eig{i} = lam_D .* a_i + b_i;            % K x nDeployments (eigenVALUES, not inverted -- see usage below)
                            end

                            %%%%%%%%%%%%%%%%% ESTIMATION %%%%%%%%%%%%%%%%%
                            total_est_start = tic;
    
                            log_msg(verbosity_level, 3, 'Begin ML estimation');
                            for strat_idx = 1:num_strats
                                strat = all_strats(strat_idx);
                                log_msg(verbosity_level, 4, 'Running ML estimation for %s | obj = %s', strat, obj_func);
        
                                ai = 1;
                                bi = 0;

                                y = pagetranspose(out_cws + scaled_n);
                                if attacker_enabled
                                    y_no_attack = pagetranspose(out_cws_no_attack + scaled_n);   % SAME noise draw, only the attacked sensor's content differs
                                end

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
    
                                    if contains(experiment, "disjoint")
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
                                    g_sq = reshape(g_tilde, S, num_antennas, nDeployments);   % S x num_antennas x nDeployments

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

                                    %% --- Subspace-projection residual test: whitened
                                    % out-of-Phi energy per sensor. T_i ~ chi^2_{K-d_sub}
                                    % under "sensor i honest" -- uses NEITHER alpha_true nor
                                    % t0_true anywhere in this computation. ---
                                    % T_i_all = zeros(1, S, nTrials, nDeployments);
                                    % for i = 1:S
                                    %     u_i = reshape(hat_u_all(:,i,:,:), K, nTrials, nDeployments);
                                    % 
                                    %     % In-Phi part (what a matching honest shape explains)
                                    %     % removed; Psi-coordinates of whatever's left.
                                    %     proj_i = pagemtimes(Phi, pagemtimes(Phi.', u_i));
                                    %     R_i_psi = pagemtimes(Psi.', u_i - proj_i);   % (K-d_sub) x nTrials x nDeployments
                                    % 
                                    %     R_i_row = reshape(R_i_psi, 1, K-d_sub, nTrials, nDeployments);
                                    %     R_i_col = reshape(R_i_psi, K-d_sub, 1, nTrials, nDeployments);
                                    %     Ki_inv_bcast = reshape(Ki_inv_all{i}, K-d_sub, K-d_sub, 1, nDeployments);
                                    % 
                                    %     % Whitened quadratic form: T_i = R_i^T * Ki^-1 * R_i
                                    %     T_i_all(1,i,:,:) = pagemtimes(pagemtimes(R_i_row, Ki_inv_bcast), R_i_col);
                                    % end

                                    T_i_all = zeros(1, S, nTrials, nDeployments);
                                    for i = 1:S
                                        u_i = reshape(hat_u_all(:,i,:,:), K, nTrials, nDeployments);

                                        proj_i = pagemtimes(Phi, pagemtimes(Phi.', u_i));
                                        R_i_psi = pagemtimes(Psi.', u_i - proj_i);   % (K-d_sub) x nTrials x nDeployments

                                        % Rotate into Rss_Psi's eigenbasis -- whitening is now an
                                        % elementwise reweight, not a full (K-d_sub)x(K-d_sub) form.
                                        R_i_v = pagemtimes(V_Rss.', R_i_psi);   % (K-d_sub) x nTrials x nDeployments
                                        weights = reshape(Ki_diag_all{i}, K-d_sub, 1, nDeployments);

                                        T_i_all(1,i,:,:) = reshape(sum((R_i_v.^2) .* weights, 1), 1,1,nTrials,nDeployments);
                                    end

                                    % Fixed, precomputed threshold -- same value for every
                                    % sensor, every trial, every deployment (Statistics and
                                    % Machine Learning Toolbox required for chi2inv).
                                    delta_fa = 1e-6;
                                    T_threshold = chi2inv(1-delta_fa, K-d_sub);
                                    flagged = T_i_all > T_threshold;   % 1 x S x nTrials x nDeployments logical

                                    %% --- Practical hat_t0 for the sign-flip check: robust
                                    % version. (1) EXCLUDES sensors already flagged by T_i --
                                    % an attacker whose shape already failed T_i should not
                                    % pollute this timing estimate. (2) Uses the MEDIAN of
                                    % independent PER-SENSOR delay estimates, not a raw sum --
                                    % exploits |B| < S/2 directly, so a minority of attackers
                                    % (even unflagged ones) cannot bias it. Own diagnostic use
                                    % only -- not the final reported t0_hat. Same circularity
                                    % caveat as T_i^align: still uses raw, unscreened DATA for
                                    % T_i itself, just no longer double-counts already-flagged
                                    % sensors on top of that.
                                    corr_per_sensor = dt * pagemtimes(D_template.', reshape(hat_u_all, K, S*nTrials*nDeployments));   % K x (S*nTrials*nDeployments)
                                    [~, best_col_per_sensor] = max(corr_per_sensor(offset_idx:end,:), [], 1);

                                    best_col_per_sensor = best_col_per_sensor + offset_idx - 1;
                                    t0_col_idx_per_sensor = reshape(best_col_per_sensor, S, nTrials, nDeployments) - 1;   % per-sensor delay guess

                                    t0_col_idx_per_trial = zeros(1, nTrials, nDeployments);
                                    for d_idx = 1:nDeployments
                                        for tr = 1:nTrials
                                            already_flagged = find(flagged(1,:,tr,d_idx));
                                            candidates = setdiff(1:S, already_flagged);
                                            if isempty(candidates)
                                                candidates = 1:S;   % fallback: everyone flagged, use all anyway
                                            end
                                            t0_col_idx_per_trial(1,tr,d_idx) = round(median(t0_col_idx_per_sensor(candidates,tr,d_idx)));
                                        end
                                    end

                                    %% --- Sign-flip defense ---
                                    % Correlate each sensor's recovered signal against the
                                    % KNOWN template rho(t-t0_true) -- not against t0_true
                                    % itself, but the one column of D_template nearest it.
                                    % Design note: using t0_true (rather than the array's own
                                    % data-derived t0_hat) sidesteps the circularity concern
                                    % discussed earlier -- a real deployment would substitute
                                    % a pre-established/robustified consensus t0 here instead.
                                    % Sign is robust to small t0 error, so this is a reasonable
                                    % simplification for this first pass. An honest sensor's
                                    % correlation is unconditionally positive; a sign-flipped
                                    % attacker's is unconditionally negative -- no calibration,
                                    % no threshold tuning needed, unlike T_i.

                                    hat_c_all = zeros(1, S, nTrials, nDeployments);
                                    for i = 1:S
                                        u_i = reshape(hat_u_all(:,i,:,:), K, nTrials, nDeployments);
                                        for d_idx = 1:nDeployments
                                            for tr = 1:nTrials
                                                rho_ref = D_template(:, t0_col_idx_per_trial(1,tr,d_idx));
                                                hat_c_all(1,i,tr,d_idx) = dt * (rho_ref.' * u_i(:,tr,d_idx));
                                            end
                                        end
                                    end

                                    sign_flagged = hat_c_all < 0;   % 1 x S x nTrials x nDeployments logical
                                    flagged = flagged | sign_flagged;   % combine with the existing T_i-based flag

                                    %% --- Coherent-replay defense (T_i^align) ---
                                    % Tests each sensor against the SINGLE consensus direction
                                    % rho(.-t0_hat) (reusing t0_col_idx_per_trial from the
                                    % sign-flip check above), not the full d_sub-dim family Phi.
                                    % A wrong-tau replay passes T_i (fools the family test) but
                                    % should FAIL this -- low correlation with the array's own
                                    % consensus timing. GLS decomposition in D_template's own
                                    % eigenbasis (U_D): E_total = whitened energy of u_i,
                                    % S_captured = whitened energy explained by the single
                                    % direction rho(.-t0_hat), T_align = E_total - S_captured.
                                    T_align_all = zeros(1, S, nTrials, nDeployments);
                                    for i = 1:S
                                        u_i = reshape(hat_u_all(:,i,:,:), K, nTrials, nDeployments);
                                        u_i_v = pagemtimes(U_D.', u_i);   % K x nTrials x nDeployments, in D_template's eigenbasis

                                        for d_idx = 1:nDeployments
                                            eig_i_d = Ki_full_eig{i}(:,d_idx);   % K x 1

                                            for tr = 1:nTrials
                                                rho_ref = D_template(:, t0_col_idx_per_trial(1,tr,d_idx));
                                                rho_ref_v = U_D.' * rho_ref;   % K x 1, same eigenbasis

                                                G_scalar = sum(rho_ref_v.^2 ./ eig_i_d);
                                                c_scalar = sum(rho_ref_v .* u_i_v(:,tr,d_idx) ./ eig_i_d);
                                                E_total  = sum(u_i_v(:,tr,d_idx).^2 ./ eig_i_d);

                                                S_captured = c_scalar^2 / G_scalar;
                                                T_align_all(1,i,tr,d_idx) = E_total - S_captured;
                                            end
                                        end
                                    end

                                    delta_fa_align = 1e-6;
                                    T_align_threshold = chi2inv(1-delta_fa_align, K-1);   % K-1 DOF: one direction removed
                                    align_flagged = T_align_all > T_align_threshold;
                                    flagged = flagged | align_flagged;

                                    %% --- Cross-sensor correlation consistency (Flagging IV) ---
                                    % Honest sensors are all copies of the SAME rho(t-bar_t0), so
                                    % they should correlate strongly with EACH OTHER, independent
                                    % of any noise-covariance model K_i. Fails differently than
                                    % T_i^align -- no whitening model needed at all -- so this is a
                                    % genuine backstop, not a repeat of the same mechanism.
                                    gram = dt * pagemtimes(pagetranspose(hat_u_all), hat_u_all);   % S x S x nTrials x nDeployments

                                    norm_sq = zeros(1,S,nTrials,nDeployments);
                                    for s = 1:S
                                        norm_sq(1,s,:,:) = gram(s,s,:,:);
                                    end
                                    norm_col = reshape(norm_sq, S,1,nTrials,nDeployments);
                                    norm_row = reshape(norm_sq, 1,S,nTrials,nDeployments);
                                    denom_corr = sqrt(pagemtimes(norm_col, norm_row));   % S x S x nTrials x nDeployments
                                    corr_mat = gram ./ denom_corr;                        % normalized pairwise correlation, in [-1,1]

                                    med_corr = zeros(1,S,nTrials,nDeployments);
                                    for i = 1:S
                                        others = setdiff(1:S,i);
                                        med_corr(1,i,:,:) = median(corr_mat(i,others,:,:), 2);
                                    end

                                    corr_delta_k = 3;   % MAD multiplier, same convention as elsewhere
                                    med_of_meds = median(med_corr, 2);
                                    mad_corr = mad(med_corr, 1, 2);

                                    corr_flagged = med_corr < (med_of_meds - corr_delta_k .* mad_corr);
                                    % flagged = flagged | corr_flagged;

                                    %% --- Individual (per-sensor) t0 and alpha estimation ---
                                    % Each sensor's isolated signal hat_u_i(t) is treated as its
                                    % OWN single-sensor estimation problem: t0_i_hat, alpha_i_hat
                                    % are computed using ONLY that sensor's own recovered data,
                                    % its own known m_i, and its own known noise model (Ki_full_eig)
                                    % -- no shared array quantity (no hat_alpha, no hat_t0) enters
                                    % anywhere. This is the "Capability A" defense: an attacker
                                    % that only falsifies transmitted content (not its own
                                    % calibration record) produces an alpha_i_hat that deviates
                                    % from the honest cluster; comparing {alpha_i_hat} across
                                    % sensors requires no external reference.
                                    %
                                    % t0_i_hat REUSES t0_col_idx_per_sensor (already computed
                                    % above via unweighted correlation search) -- the ML argmax
                                    % location over tau is unaffected by GLS whitening, since
                                    % m_i^2 > 0 scales the design vector uniformly across all
                                    % candidate tau and does not shift the argmax.
                                    alpha_i_hat = zeros(1, S, nTrials, nDeployments);
                                    t0_i_hat    = zeros(1, S, nTrials, nDeployments);
                                    G_scalar_all = zeros(1, S, nTrials, nDeployments);

                                    for i = 1:S
                                        u_i = reshape(hat_u_all(:,i,:,:), K, nTrials, nDeployments);
                                        u_i_v = pagemtimes(U_D.', u_i);   % K x nTrials x nDeployments, in D_template's eigenbasis

                                        for d_idx = 1:nDeployments
                                            eig_i_d = Ki_full_eig{i}(:,d_idx);   % K x 1
                                            mi_val = mi_5d(1,i,1,1,d_idx);

                                            for tr = 1:nTrials
                                                col_idx = t0_col_idx_per_sensor(i,tr,d_idx);
                                                t0_i_hat(1,i,tr,d_idx) = t(col_idx);

                                                rho_ref = D_template(:, col_idx);
                                                rho_ref_v = U_D.' * rho_ref;   % K x 1, same eigenbasis

                                                % GLS fit of u_i ~ (alpha_i * m_i^2) * rho_ref + noise,
                                                % whitened by this sensor's OWN known noise model.
                                                G_scalar = sum(rho_ref_v.^2 ./ eig_i_d);
                                                c_scalar = sum(rho_ref_v .* u_i_v(:,tr,d_idx) ./ eig_i_d);

                                                alpha_i_hat(1,i,tr,d_idx) = c_scalar / (G_scalar * mi_val^2);
                                                G_scalar_all(1,i,tr,d_idx) = G_scalar * mi_val^4;   % precision of alpha_i_hat itself (chain rule through the mi^2 scaling)
                                            end
                                        end
                                    end

                                    k_trim = 1;
                                    sorted_alpha = sort(alpha_i_hat, 2);
                                    alpha_estimates_trimmed = mean(sorted_alpha(1,k_trim+1:end-k_trim,:,:), 2);
                                    t0_estimates_trimmed = reshape(t(t0_col_idx_per_trial), 1,1,nTrials,nDeployments);

                                    %% --- Undefended array-level presence check (Q statistic) ---
                                    % Cheap, no per-sensor flagging, no |B|<S/2 assumption needed --
                                    % just asks "does this array's spread of independent alpha_i_hat
                                    % look larger than pure noise should produce." bar_alpha cancels
                                    % out (weighted variance around the array's OWN pooled mean), so
                                    % this has a clean, exact, absolute chi^2_{S-1} threshold -- the
                                    % legitimate version of the Lambda(empty-set) gate idea from
                                    % earlier, which failed only because it lacked this self-
                                    % consistency property.
                                    G_pooled_num = sum(G_scalar_all .* alpha_i_hat, 2);
                                    G_pooled_den = sum(G_scalar_all, 2);
                                    alpha_pooled = G_pooled_num ./ G_pooled_den;   % 1 x 1 x nTrials x nDeployments

                                    Q_stat = sum(G_scalar_all .* (alpha_i_hat - alpha_pooled).^2, 2);   % 1 x 1 x nTrials x nDeployments
                                    delta_fa_Q = 1e-3;
                                    Q_threshold = chi2inv(1-delta_fa_Q, S-1);

                                    array_suspected = Q_stat > Q_threshold;   % 1 x 1 x nTrials x nDeployments logical

                                    %% --- Sensor-vs-sensor amplitude consistency (Capability A defense) ---
                                    % Honest sensors' independent alpha_i_hat estimates should all
                                    % cluster near the SAME true bar_alpha -- an attacker that only
                                    % falsifies transmitted content (m_i on file stays correct) will
                                    % show alpha_i_hat displaced by exactly its amplitude lie factor.
                                    % Relative (median/MAD) test, same reasoning as every other
                                    % robust check in this pipeline: needs |B| < S/2, no external
                                    % reference or shared array quantity required.
                                    alpha_delta_k = 10;   % MAD multiplier, same convention as elsewhere
                                    med_alpha = median(alpha_i_hat, 2);   % 1 x 1 x nTrials x nDeployments
                                    mad_alpha = mad(alpha_i_hat, 1, 2);   % median absolute deviation across sensors

                                    amplitude_flagged = abs(alpha_i_hat - med_alpha) > (alpha_delta_k .* mad_alpha);
                                    % flagged = flagged | amplitude_flagged;

                                    honest_idx = setdiff(1:S, attacker_idx);

                                    %% --- DIAGNOSTIC: isolate Nw-term and N0-term scaling in K_i ---
                                    % diag_sensor = setdiff(1:S, attacker_idx); diag_sensor = diag_sensor(1);   % an honest sensor
                                    % diag_dep = 1;
                                    % 
                                    % Pi_bcast      = reshape(Pi_all{diag_sensor}, r_dim, Mtot, 1, nDeployments);
                                    % breve_bcast   = reshape(breve_g_all{diag_sensor}, 1, r_dim, 1, nDeployments);
                                    % norm_sq_bcast = reshape(sum(breve_g_all{diag_sensor}.^2,1), 1,1,1,nDeployments);
                                    % 
                                    % % --- Test A: sensor-noise-only (zero antenna noise) ---
                                    % [~, ui_noise_only] = mf_integral_fft(scaled_w(:,1:S,:,:,:), mi_5d .* sensor_signal(t, Tp, norm_fact), 1, 1, K, dt, Tp);
                                    % y_A = pagetranspose(sum(g_tilde .* ui_noise_only, 2));           % NO scaled_n added
                                    % y_sq_A = reshape(y_A, K, num_antennas, nTrials, nDeployments);
                                    % y_R_A = permute(cat(2, real(y_sq_A), imag(y_sq_A)), [2 1 3 4]);
                                    % 
                                    % hat_u_A = pagemtimes(breve_bcast, pagemtimes(Pi_bcast, y_R_A)) ./ norm_sq_bcast;
                                    % hat_u_A = reshape(hat_u_A(:,:,:,diag_dep), K, nTrials);
                                    % Ri_psi_A = Psi.' * (hat_u_A - Phi*(Phi.'*hat_u_A));
                                    % 
                                    % mi_val = mi_5d(1,diag_sensor,1,1,diag_dep);
                                    % ratio_A = trace(Ri_psi_A*Ri_psi_A.'/nTrials) / trace(mi_val^2 * Nw_over_2 * Rss_Psi);
                                    % fprintf('Nw term: empirical/predicted = %.4f   (1/dt = %.4f)\n', ratio_A, 1/dt);
                                    % 
                                    % % --- Test B: antenna-noise-only (zero sensor noise / signal) ---
                                    % y_B = pagetranspose(scaled_n);
                                    % y_sq_B = reshape(y_B, K, num_antennas, nTrials, nDeployments);
                                    % y_R_B = permute(cat(2, real(y_sq_B), imag(y_sq_B)), [2 1 3 4]);
                                    % 
                                    % hat_u_B = pagemtimes(breve_bcast, pagemtimes(Pi_bcast, y_R_B)) ./ norm_sq_bcast;
                                    % hat_u_B = reshape(hat_u_B(:,:,:,diag_dep), K, nTrials);
                                    % Ri_psi_B = Psi.' * (hat_u_B - Phi*(Phi.'*hat_u_B));
                                    % 
                                    % norm_sq_i = sum(breve_g_all{diag_sensor}(:,diag_dep).^2);
                                    % ratio_B = trace(Ri_psi_B*Ri_psi_B.'/nTrials) / trace((N0_over_2/norm_sq_i)*eye(K-d_sub));
                                    % fprintf('N0 term: empirical/predicted = %.4f   (1/dt = %.4f)\n', ratio_B, 1/dt);
                                    
                                    G_R = cat(2, real(g_sq), imag(g_sq));
                                    G_R = permute(G_R, [2,1,3]);

                                    W_bcast = reshape(W, Mtot, Mtot, 1, nDeployments);  % broadcast over nTrials

                                    z = pagemtimes(W_bcast, y_R);   % 2M x K x nTrials x nDeployments, decorrelated signal+noise
                                    z_bs = pagemtimes(W_bcast, y_R_bs);

                                    mu = m.^2;
                                    WG_Rmu = reshape(pagemtimes(W, pagemtimes(G_R, reshape(mu, S, 1, nDeployments))), Mtot, nDeployments);  % Mtot x nDeployments

                                    mf_with_z_sum = zeros(1, 1, K, nTrials, nDeployments);   % running z_chi(t0)
                                    num   = zeros(1,1,nTrials,nDeployments);
                                    denom = zeros(1,1,nTrials,nDeployments);
                                    Qn_matrix_all = cell(Mtot, 1);   % cache for reuse in alpha estimation

                                    for m_idx = 1:Mtot
                                        lambda_m = reshape(lambda_vals(m_idx,:), 1, 1, 1, nDeployments);
                                        mag_sqr_H_m = Hm_arr(omega, lambda_m);

                                        [~, Qn_matrix_m] = get_time_domain(mag_sqr_H_m,dt,nDeployments);
                                        resh_Qn_m_5d = reshape(Qn_matrix_m,K,K,1,1,nDeployments);
                                        Qn_matrix_all{m_idx} = Qn_matrix_m;   % <-- cache it

                                        b_m = reshape(WG_Rmu(m_idx,:), 1, 1, 1, nDeployments);
                                        Omega_m_t0_search = reshape(sum(b_m .* R00_t0_search,2), K, 1, K, 1, nDeployments);

                                        % Branch-m decorrelated observation z_m(t)
                                        z_m = z(m_idx,:,:,:);
                                        resh_z_m_5d = reshape(z_m, 1, K, 1, nTrials, nDeployments);

                                        % Accumulate z_chi(t0) and Gamma_chi(t0) contributions from branch m
                                        mf_with_z_sum = mf_with_z_sum + dt*dt*pagemtimes(pagemtimes(resh_z_m_5d, resh_Qn_m_5d), Omega_m_t0_search);
                                    end

                                    mf_with_z_sum_bs = zeros(1,1,K,nTrials,nDeployments);
                                    for m_idx = 1:Mtot
                                        b_m = reshape(WG_Rmu(m_idx,:),1,1,1,nDeployments);
                                        Omega_m_t0_search = reshape(sum(b_m .* R00_t0_search,2), K,1,K,1,nDeployments);
                                        z_m_bs = reshape(z_bs(m_idx,:,:,:),1,K,1,nTrials,nDeployments);
                                        mf_with_z_sum_bs = mf_with_z_sum_bs + dt*dt*pagemtimes(pagemtimes(z_m_bs, reshape(Qn_matrix_all{m_idx},K,K,1,1,nDeployments)), Omega_m_t0_search);
                                    end
                                    [~,I_bs] = max(mf_with_z_sum_bs(:,:,offset_idx:end,:,:),[],3);
                                    I_bs = I_bs + offset_idx - 1;
                                    t0_estimates_bs = reshape((I_bs-1)*dt,1,1,nTrials,nDeployments);

                                    % Get time domain matrix for Qn.
                                    [max_y_vals,I] = max(mf_with_z_sum(:,:,offset_idx:end,:,:), [], 3);
                                    I = I + offset_idx - 1;
                                    t0_estimates_for_plot = reshape((I-1)*dt, 1, 1, nTrials, nDeployments);
                                    t0_estimates_for_alpha = t0_true*ones(1,1,nTrials,nDeployments);

                                    % Recompute the RAW (unweighted) pulse correlation at t0_estimates_for_alpha.
                                    % No mi weighting here -- mu = mi^2 is already folded into b_m via W*G_R*mu.
                                    Af = fft(sensor_signal(t-t0_estimates_for_alpha, Tp, norm_fact), N, 1);
                                    lag0 = round(Tp/dt);
                                    Y_tensor = dt*ifft(Af .* mfTemplateFFT_raw, [], 1);
                                    Rss_tensor = Y_tensor(lag0+1 : (lag0+K), :, :, :);   % K x 1 x nTrials x nDeployments

                                    for m_idx = 1:Mtot
                                        Qn_matrix_m = Qn_matrix_all{m_idx};
                                        resh_Qn_m = reshape(Qn_matrix_m, K, K, 1, nDeployments);

                                        % Branch-m correlator kernel at the estimated t0: Omega_m(t0_hat) = b_m * Rss(t, t0_hat)
                                        b_m = reshape(WG_Rmu(m_idx,:), 1, 1, 1, nDeployments);
                                        Omega_m = sum(b_m .* Rss_tensor, 2);
                                        resh_Omega_m = reshape(Omega_m, 1, K, nTrials, nDeployments);

                                        % Branch-m observation
                                        z_m = z(m_idx,:,:,:);

                                        num_m   = dt*dt*pagemtimes(pagemtimes(z_m, resh_Qn_m), pagetranspose(resh_Omega_m));
                                        denom_m = dt*dt*pagemtimes(pagemtimes(resh_Omega_m, resh_Qn_m), pagetranspose(resh_Omega_m));

                                        num   = num   + num_m;
                                        denom = denom + denom_m;
                                    end
                                    alpha_estimates = num ./ denom;

                                    num_bs = zeros(1,1,nTrials,nDeployments); denom_bs = num_bs;
                                    for m_idx = 1:Mtot
                                        b_m = reshape(WG_Rmu(m_idx,:),1,1,1,nDeployments);
                                        resh_Omega_bs = reshape(b_m .* Rss_tensor,1,K,nTrials,nDeployments);
                                        z_m_bs = z_bs(m_idx,:,:,:);
                                        Qn_bc = reshape(Qn_matrix_all{m_idx},K,K,1,nDeployments);
                                        num_bs = num_bs + dt*dt*pagemtimes(pagemtimes(z_m_bs,Qn_bc),pagetranspose(resh_Omega_bs));
                                        denom_bs = denom_bs + dt*dt*pagemtimes(pagemtimes(resh_Omega_bs,Qn_bc),pagetranspose(resh_Omega_bs));
                                    end
                                    alpha_estimates_bs = num_bs ./ denom_bs;

                                    % --- Baseline: full S-sensor array, as if the attacker never
                                    % hijacked anyone. Reuses build_null_geometry/estimate_with_
                                    % geometry with an EMPTY exclusion set -- A=[] means every
                                    % sensor is retained, so this reduces to the same eigen-
                                    % decomposition the ordinary (undefended) pipeline already
                                    % runs, just fed the pristine y_no_attack instead of y.
                                    alpha_estimates_baseline = zeros(size(alpha_estimates));
                                    t0_estimates_baseline    = zeros(size(t0_estimates_for_plot));
                                    for d_idx = 1:nDeployments
                                        geom_baseline = build_null_geometry([], mi_5d, g_tilde, gamma_w, gamma_n, Hm_arr, omega, dt, K, S, d_idx);
                                        [alpha_b, t0_b, ~] = estimate_with_geometry(geom_baseline, y_no_attack, K, N, Tp, ...
                                            norm_fact, t, t0_true, mfTemplateFFT_raw, D_template, 1:nTrials, d_idx, dt);
                                        alpha_estimates_baseline(1,1,:,d_idx) = alpha_b;
                                        t0_estimates_baseline(1,1,:,d_idx) = t0_b;
                                    end

                                    % Lambda(empty-set): captured coherent energy of the
                                    % UNDEFENDED, full-array estimate -- free, already have
                                    % both ingredients. This is the baseline every candidate
                                    % nulling gets compared against.
                                    Lambda_empty = (alpha_estimates.^2) .* denom;

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

                                        % --- Oracle ceiling: null the TRUE attacker_idx directly,
                                        % no detection, no validation -- upper bound on what
                                        % identification+nulling could ever achieve. Reuses the
                                        % exact same geometry/estimation machinery as the real
                                        % defense, just fed the ground-truth exclusion set.
                                        alpha_estimates_oracle = zeros(size(alpha_estimates));
                                        t0_estimates_oracle    = zeros(size(t0_estimates_for_plot));
                                        for d_idx = 1:nDeployments
                                            geom_oracle = build_null_geometry(attacker_idx, mi_5d, g_tilde, gamma_w, gamma_n, Hm_arr, omega, dt, K, S, d_idx);
                                            [alpha_o, t0_o, ~] = estimate_with_geometry(geom_oracle, y, K, N, Tp, ...
                                                norm_fact, t, t0_true, mfTemplateFFT_raw, D_template, 1:nTrials, d_idx, dt);
                                            alpha_estimates_oracle(1,1,:,d_idx) = alpha_o;
                                            t0_estimates_oracle(1,1,:,d_idx) = t0_o;
                                        end
                                    end

                                    alpha_final = alpha_estimates;
                                    t0_final = t0_estimates_for_plot;

                                    if use_greedy_validation
                                        % --- Greedy, Lambda-validated: genuinely per-trial, since
                                        % accept/reject depends on each trial's own noise realization
                                        % (see prior discussion -- grouping here would be invalid). ---
                                        geom_cache = containers.Map('KeyType','char','ValueType','any');

                                        accepted_count = 0;
                                        rejected_count = 0;
                                        tp_total = 0;
                                        fp_total = 0;
                                        fn_total = 0;

                                        for d_idx = 1:nDeployments
                                            for tr = 1:nTrials
                                                F = excluded_sets{tr,d_idx};
                                                if isempty(F)
                                                    continue   % nothing flagged -- keep undefended estimate
                                                end

                                                % Order THIS trial's candidates by ITS OWN T_i --
                                                % fully per-trial, no averaging across trials.
                                                amp_deviation = abs(alpha_i_hat(1,F,tr,d_idx) - med_alpha(1,1,tr,d_idx)) ./ mad_alpha(1,1,tr,d_idx);
                                                corr_deviation = (med_of_meds(1,1,tr,d_idx) - med_corr(1,F,tr,d_idx)) ./ mad_corr(1,1,tr,d_idx);
                                                Ti_this = T_i_all(1,F,tr,d_idx) / T_threshold ...
                                                        + T_align_all(1,F,tr,d_idx) / T_align_threshold ...
                                                        + amp_deviation ...
                                                        + max(corr_deviation, 0) ...
                                                        + 50 * array_suspected(1,1,tr,d_idx) ...   % Q: corroborating group-level signal, not a standalone detector
                                                        + 1e6 * sign_flagged(1,F,tr,d_idx);
                                                [~, order] = sort(Ti_this(:), 'descend');
                                                F_sorted = F(order.');

                                                current_A      = [];
                                                current_alpha  = alpha_estimates(1,1,tr,d_idx);
                                                current_t0     = t0_estimates_for_plot(1,1,tr,d_idx);
                                                current_Lambda = Lambda_empty(1,1,tr,d_idx);

                                                for c = F_sorted
                                                    try_A = sort([current_A, c]);   % sorted -> order-independent cache key
                                                    key = sprintf('%d_%s', d_idx, mat2str(try_A));

                                                    if isKey(geom_cache, key)
                                                        geom = geom_cache(key);
                                                    else
                                                        geom = build_null_geometry(try_A, mi_5d, g_tilde, gamma_w, gamma_n, Hm_arr, omega, dt, K, S, d_idx);
                                                        geom_cache(key) = geom;
                                                    end

                                                    [try_alpha, try_t0, try_Lambda] = estimate_with_geometry(geom, y, K, N, Tp, ...
                                                        norm_fact, t, t0_true, mfTemplateFFT_raw, D_template, tr, d_idx, dt);

                                                    if try_Lambda > current_Lambda
                                                        current_A      = try_A;
                                                        current_alpha  = try_alpha;
                                                        current_t0     = try_t0;
                                                        current_Lambda = try_Lambda;
                                                        accepted_count = accepted_count + 1;
                                                    else
                                                        rejected_count = rejected_count + 1;
                                                    end
                                                end

                                                true_positive  = numel(intersect(current_A, attacker_idx));
                                                false_positive = numel(setdiff(current_A, attacker_idx));
                                                false_negative = numel(setdiff(attacker_idx, current_A));
                                                tp_total = tp_total + true_positive;
                                                fp_total = fp_total + false_positive;
                                                fn_total = fn_total + false_negative;

                                                if ~isempty(current_A)
                                                    alpha_final(1,1,tr,d_idx) = current_alpha;
                                                    t0_final(1,1,tr,d_idx) = current_t0;
                                                end
                                            end
                                        end

                                        log_msg(verbosity_level, 4, 'Greedy validation: %d accepted, %d rejected (%d distinct null-geometries built and cached)', ...
                                            accepted_count, rejected_count, geom_cache.Count);
                                        log_msg(verbosity_level, 4, 'Detection accuracy: TP=%d, FP=%d, FN=%d', tp_total, fp_total, fn_total);
                                    else
                                        % --- Face value: no per-trial decision to make, so group
                                        % trials by their EXACT flagged set and process each group
                                        % in ONE batched call -- both geometry-building AND
                                        % estimation are shared across every trial in the group. ---
                                        num_geometries_built = 0;
                                        tp_total = 0;
                                        fp_total = 0;
                                        fn_total = 0;

                                        for d_idx = 1:nDeployments
                                            [unique_sets, group_idx] = unique_cell_sets(excluded_sets(:,d_idx));

                                            for g = 1:numel(unique_sets)
                                                A = unique_sets{g};
                                                if isempty(A)
                                                    continue   % nothing flagged -- keep undefended estimate
                                                end

                                                trial_members = find(group_idx == g);

                                                geom = build_null_geometry(A, mi_5d, g_tilde, gamma_w, gamma_n, Hm_arr, omega, dt, K, S, d_idx);
                                                num_geometries_built = num_geometries_built + 1;

                                                [alpha_A, t0_A, ~] = estimate_with_geometry(geom, y, K, N, Tp, ...
                                                    norm_fact, t, t0_true, mfTemplateFFT_raw, D_template, trial_members, d_idx, dt);

                                                true_positive  = numel(intersect(A, attacker_idx));
                                                false_positive = numel(setdiff(A, attacker_idx));
                                                false_negative = numel(setdiff(attacker_idx, A));
                                                tp_total = tp_total + true_positive * numel(trial_members);
                                                fp_total = fp_total + false_positive * numel(trial_members);
                                                fn_total = fn_total + false_negative * numel(trial_members);

                                                alpha_final(1,1,trial_members,d_idx) = alpha_A;
                                                t0_final(1,1,trial_members,d_idx) = t0_A;
                                            end
                                        end

                                        log_msg(verbosity_level, 4, 'Face-value nulling: %d distinct null-geometries built across all deployments', num_geometries_built);
                                        log_msg(verbosity_level, 4, 'Detection accuracy: TP=%d, FP=%d, FN=%d', tp_total, fp_total, fn_total);

                                        honest_idx_diag = setdiff(1:S, attacker_idx);
                                        fp_Ti     = mean(T_i_all(1,honest_idx_diag,:,:) > T_threshold, 'all');
                                        fp_sign   = mean(sign_flagged(1,honest_idx_diag,:,:), 'all');
                                        fp_align  = mean(align_flagged(1,honest_idx_diag,:,:), 'all');
                                        fp_amplitude = mean(amplitude_flagged(1,honest_idx_diag,:,:), 'all');
                                        fp_corr = mean(corr_flagged(1,honest_idx_diag,:,:), 'all');
                                        fp_combined = mean(flagged(1,honest_idx_diag,:,:), 'all');

                                        log_msg(verbosity_level, 4, 'FP rate -- T_i: %.4f | sign: %.4f | align: %.4f | amplitude: %0.4f | corr: %.4f | combined: %.4f', ...
                                            fp_Ti, fp_sign, fp_align, fp_amplitude, fp_corr, fp_combined);
                                    end

                                    alpha_estimates = alpha_final;
                                    t0_estimates_for_plot = t0_final;
                                end

                                % Compute empirical variance.
                                rho_empirical_var(strat_idx,agent_db_idx,channel_db_idx,1,scheme_idx,pivot_idx,save_dim,:) = var(alpha_estimates,0,3);
                                rho_empirical_var(strat_idx,agent_db_idx,channel_db_idx,2,scheme_idx,pivot_idx,save_dim,:) = var(t0_estimates_for_plot,0,3);

                                % Compute empirical mse.
                                rho_empirical_mse(strat_idx,agent_db_idx,channel_db_idx,1,scheme_idx,pivot_idx,save_dim,:) = mean((alpha_estimates - alpha_true).^2,3);
                                rho_empirical_mse(strat_idx,agent_db_idx,channel_db_idx,2,scheme_idx,pivot_idx,save_dim,:) = mean((t0_estimates_for_plot - t0_true).^2,3);

                                % Compute UNDEFENDED empirical mse -- only under attack, per request.
                                if attacker_enabled
                                    rho_empirical_mse_baseline(strat_idx,agent_db_idx,channel_db_idx,1,scheme_idx,pivot_idx,save_dim,:) = mean((alpha_estimates_baseline - alpha_true).^2,3);
                                    rho_empirical_mse_baseline(strat_idx,agent_db_idx,channel_db_idx,2,scheme_idx,pivot_idx,save_dim,:) = mean((t0_estimates_baseline - t0_true).^2,3);

                                    rho_empirical_mse_undefended(strat_idx,agent_db_idx,channel_db_idx,1,scheme_idx,pivot_idx,save_dim,:) = mean((alpha_estimates_undefended - alpha_true).^2,3);
                                    rho_empirical_mse_undefended(strat_idx,agent_db_idx,channel_db_idx,2,scheme_idx,pivot_idx,save_dim,:) = mean((t0_estimates_undefended - t0_true).^2,3);
                                
                                    rho_empirical_mse_oracle(strat_idx,agent_db_idx,channel_db_idx,1,scheme_idx,pivot_idx,save_dim,:) = mean((alpha_estimates_oracle - alpha_true).^2,3);
                                    rho_empirical_mse_oracle(strat_idx,agent_db_idx,channel_db_idx,2,scheme_idx,pivot_idx,save_dim,:) = mean((t0_estimates_oracle - t0_true).^2,3);

                                    rho_empirical_mse_bs(strat_idx,agent_db_idx,channel_db_idx,1,scheme_idx,pivot_idx,save_dim,:) = mean((alpha_estimates_bs - alpha_true).^2,3);
                                    rho_empirical_mse_bs(strat_idx,agent_db_idx,channel_db_idx,2,scheme_idx,pivot_idx,save_dim,:) = mean((t0_estimates_bs - t0_true).^2,3);

                                    rho_empirical_mse_trimmed(strat_idx,agent_db_idx,channel_db_idx,1,scheme_idx,pivot_idx,save_dim,:) = mean((alpha_estimates_trimmed - alpha_true).^2,3);
                                    rho_empirical_mse_trimmed(strat_idx,agent_db_idx,channel_db_idx,2,scheme_idx,pivot_idx,save_dim,:) = mean((t0_estimates_trimmed - t0_true).^2,3);
                                end

                                % Compute empirical bias.
                                % rho_empirical_bias(strat_idx,agent_db_idx,channel_db_idx,1,scheme_idx,pivot_idx,save_dim,:) = mean(alpha_estimates,3);
                                % rho_empirical_bias(strat_idx,agent_db_idx,channel_db_idx,2,scheme_idx,pivot_idx,save_dim,:) = mean(t0_estimates_for_plot,3);

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

    params_latex = ["\hat{\alpha}","\hat{t}_0"];
    params_text = ["alpha","t0"];

    RGB = orderedcolors("gem12");
    H = compose("#%02X%02X%02X",round(RGB*255));
    colors = cellstr(H);
    colorMap = containers.Map(cellstr([unique(iter_arr)]), colors(1:numel([unique(iter_arr)])));

    exclude = all_strats(contains(all_strats, "hetero.") & contains(all_strats, "prop. scaling")); % ["None"];

    for agent_db_idx = 1:length(agent_db_values)
        selected_strat_idxs = find((1:numel(all_strats)) .* ~ismember(all_strats,exclude).');
        selected_strats = all_strats(selected_strat_idxs);
        for param_idx = 1:2
            % normalize performance metrics
            if param_idx == 1
                norm_coeff = alpha_true^2;
            else
                norm_coeff = t0_true^2;
            end

            avg_dep_crlb(selected_strat_idxs,agent_db_idx,:,param_idx,:,:) = avg_dep_crlb(selected_strat_idxs,agent_db_idx,:,param_idx,:,:) ./ norm_coeff;
            avg_dep_var(selected_strat_idxs,agent_db_idx,:,param_idx,:,:) = avg_dep_var(selected_strat_idxs,agent_db_idx,:,param_idx,:,:) ./ norm_coeff;
            avg_dep_mse(selected_strat_idxs,agent_db_idx,:,param_idx,:,:) = avg_dep_mse(selected_strat_idxs,agent_db_idx,:,param_idx,:,:) ./ norm_coeff;
            if attacker_enabled
                avg_dep_mse_undefended(selected_strat_idxs,agent_db_idx,:,param_idx,:,:) = avg_dep_mse_undefended(selected_strat_idxs,agent_db_idx,:,param_idx,:,:) ./ norm_coeff;
                avg_dep_mse_oracle(selected_strat_idxs,agent_db_idx,:,param_idx,:,:) = avg_dep_mse_oracle(selected_strat_idxs,agent_db_idx,:,param_idx,:,:) ./ norm_coeff;
                avg_dep_mse_baseline(selected_strat_idxs,agent_db_idx,:,param_idx,:,:) = avg_dep_mse_baseline(selected_strat_idxs,agent_db_idx,:,param_idx,:,:) ./ norm_coeff;
                avg_dep_mse_trimmed(selected_strat_idxs,agent_db_idx,:,param_idx,:,:) = avg_dep_mse_trimmed(selected_strat_idxs,agent_db_idx,:,param_idx,:,:) ./ norm_coeff;
                avg_dep_mse_bs(selected_strat_idxs,agent_db_idx,:,param_idx,:,:) = avg_dep_mse_bs(selected_strat_idxs,agent_db_idx,:,param_idx,:,:) ./ norm_coeff;
                avg_dep_crlb_oracle(selected_strat_idxs,agent_db_idx,:,param_idx,:,:) = avg_dep_crlb_oracle(selected_strat_idxs,agent_db_idx,:,param_idx,:,:) ./ norm_coeff;
            end

            % Collect all plot values for y-axis scaling.
            all_vals = [avg_dep_crlb(selected_strat_idxs,agent_db_idx,:,param_idx,:,1:sensor_dimension);
                        avg_dep_var(selected_strat_idxs,agent_db_idx,:,param_idx,:,1:sensor_dimension);
                        avg_dep_mse(selected_strat_idxs,agent_db_idx,:,param_idx,:,1:sensor_dimension)
                        ];
            if attacker_enabled
                all_vals = [all_vals; avg_dep_mse_undefended(selected_strat_idxs,agent_db_idx,:,param_idx,:,1:sensor_dimension);
                            avg_dep_mse_oracle(selected_strat_idxs,agent_db_idx,:,param_idx,:,1:sensor_dimension);
                            avg_dep_mse_baseline(selected_strat_idxs,agent_db_idx,:,param_idx,:,1:sensor_dimension);
                            avg_dep_mse_trimmed(selected_strat_idxs,agent_db_idx,:,param_idx,:,1:sensor_dimension);
                            avg_dep_mse_bs(selected_strat_idxs,agent_db_idx,:,param_idx,:,1:sensor_dimension);
                            avg_dep_crlb_oracle(selected_strat_idxs,agent_db_idx,:,param_idx,:,1:sensor_dimension)];
            end

            for scheme_idx = 1:length(selected_schemes)
                scheme = selected_schemes(scheme_idx);

                count_vals = sensor_vals;

                for constraint_idx = 1:length(constraints)
                    constraint = constraints(constraint_idx);
                    data_idx_seq = selected_strat_idxs(find(contains(selected_strats, constraint)).');

                    if constraint == "indiv."
                        figLgd = figure;
                        axLgd = gca;
                        plot_w = 450;
                        set(figLgd,'Position',[100,100,plot_w,plot_w])
                        hold(axLgd, 'on');

                        lgd = [];

                        for count_idx = data_idx_seq
                            strat = all_strats(count_idx);
                            split_str = split(strat);
                            iter_key = join(split_str(2:end), " ");

                            lgd = [lgd iter_key];

                            plot(axLgd, nan,nan,'color',colorMap(iter_key),'LineWidth',2);
                        end

                        plot(axLgd, nan,nan,'x','color','black');
                        plot(axLgd, nan,nan,'^','color','black');
                        plot(axLgd, nan,nan,'o','color','black');
                        legend_entries = ["MSE", "VAR", "CRLB"];
                        if attacker_enabled
                            plot(axLgd, nan,nan,'s','color','black');
                            plot(axLgd, nan,nan,'s','color','black');
                            legend_entries = [legend_entries, "MSE (undefended)", "MSE (oracle)"];
                        end

                        lgdObj = legend(axLgd, [lgd, legend_entries]);

                        lgdObj = legend(axLgd, [lgd, "MSE", "VAR", "CRLB"]);
                        lgdObj.Location = 'none';
                        lgdObj.Units = 'normalized';
                        lgdObj.Position(1) = (1 - lgdObj.Position(3)) / 2;  % center horizontally
                        lgdObj.Position(2) = (1 - lgdObj.Position(4)) / 2;  % center vertically

                        lgdObj.FontSize = 15;        % increase text size
                        lgdObj.ItemTokenSize = [30 18];  % increase marker/line size [width height]
                        axLgd.Color = 'none';
                        axLgd.XColor = 'none';
                        axLgd.YColor = 'none';
                        figLgd.Color = 'white';
                        drawnow;
                    end

                    fig1 = figure;
                    plot_w = 450;
                    set(fig1,'Position',[100,100,plot_w,plot_w])
                    ax = gca;
                    hold(ax, 'on');

                    x_axis_series = channel_db_values(scheme_idx,:,agent_db_idx);
                    plot_line_width = 1.5;

                    for count_idx = 1:length(count_vals)
                        for strat_idx = data_idx_seq
                            strat = all_strats(strat_idx);
                            split_str = split(strat);
                            iter_key = join(split_str(2:end), " ");

                            plot(ax, x_axis_series,(squeeze(avg_dep_mse(strat_idx,agent_db_idx,:,param_idx,scheme_idx,count_idx,1))),'-x','color',colorMap(iter_key),'LineWidth', plot_line_width)
                            % plot(ax, x_axis_series,(squeeze(avg_dep_var(strat_idx,agent_db_idx,:,param_idx,scheme_idx,count_idx,1))),'-^','color',colorMap(iter_key),'LineWidth', plot_line_width)
                            % plot(ax, x_axis_series,(squeeze(avg_dep_crlb(strat_idx,agent_db_idx,:,param_idx,scheme_idx,count_idx,1))),'--o','color',colorMap(iter_key),'LineWidth', plot_line_width)
                            if attacker_enabled
                                plot(ax, x_axis_series,(squeeze(avg_dep_mse_undefended(strat_idx,agent_db_idx,:,param_idx,scheme_idx,count_idx,1))),':s','color','red','LineWidth', plot_line_width)
                                % plot(ax, x_axis_series,(squeeze(avg_dep_mse_oracle(strat_idx,agent_db_idx,:,param_idx,scheme_idx,count_idx,1))),':s','color','green','LineWidth', plot_line_width)
                                plot(ax, x_axis_series,(squeeze(avg_dep_mse_baseline(strat_idx,agent_db_idx,:,param_idx,scheme_idx,count_idx,1))),':s','color','blue','LineWidth', plot_line_width)
                                plot(ax, x_axis_series,(squeeze(avg_dep_mse_trimmed(strat_idx,agent_db_idx,:,param_idx,scheme_idx,count_idx,1))),':s','color','cyan','LineWidth', plot_line_width)
                                plot(ax, x_axis_series,(squeeze(avg_dep_mse_bs(strat_idx,agent_db_idx,:,param_idx,scheme_idx,count_idx,1))),':s','color','magenta','LineWidth', plot_line_width)
                                % plot(ax, x_axis_series,(squeeze(avg_dep_crlb_oracle(strat_idx,agent_db_idx,:,param_idx,scheme_idx,count_idx,1))),':s','color','magenta','LineWidth', plot_line_width)
                            end
                        end
                    end

                    xlabel('Channel SNR (dB)')
                    ylabel(" ")

                    title('$$'+params_latex(param_idx)+'$$','Interpreter','latex')

                    ax.FontSize = 15;          % font size
                    ax.YScale = 'log';          % log scale if needed
                    grid(ax, 'on');             % turn on grid

                    y_lower = 10^(-0.25)*min(all_vals,[],"all");
                    y_upper = 1.05*max(all_vals,[],"all");
                    ylim([y_lower,y_upper]);

                    ax.XTick = channel_snr;
                end
            end
        end
    end
end

%% Functions
function log_msg(verbosity_level, level, fmt, varargin)
% level 1 = major section (scheme/agent-SNR/experiment boundaries)
% level 2 = sub-step (channel SNR point, per-block timing)
% level 3 = detail/diagnostic (per-strategy, T_i/geometry/detection stats)
% level 4 = fine-grained (per-sensor, per-candidate detail)
% level 5 = trace (innermost loop, per-trial/per-iteration detail)
if level > verbosity_level
    return
end
indent = repmat('  ', 1, level-1);
switch level
    case 1
        prefix = sprintf('\n%s=== ', indent);
        suffix = ' ===';
    case 2
        prefix = sprintf('%s-- ', indent);
        suffix = '';
    case 3
        prefix = sprintf('%s.. ', indent);
        suffix = '';
    case 4
        prefix = sprintf('%s.... ', indent);
        suffix = '';
    case 5
        prefix = sprintf('%s...... ', indent);
        suffix = '';
    otherwise
        prefix = sprintf('%s   ', indent);
        suffix = '';
end
fprintf([prefix fmt suffix '\n'], varargin{:});
end

function geom = build_null_geometry(A, mi_5d, g_tilde, gamma_w, gamma_n, Hm_arr, omega, dt, K, S, d_idx)
% Everything needed to null set A and estimate for deployment d_idx -- depends
% ONLY on A and d_idx (and gamma_w/gamma_n, fixed for the current channel-SNR
% point). NEVER depends on trial data -- safe to cache and reuse across every
% trial and every greedy step that happens to test this same A.

num_antennas = size(g_tilde,3);
Mtot_full = 2*num_antennas;
retained = setdiff(1:S, A);

if numel(A) >= Mtot_full
    error('Cannot null %d sensors with only %d real antenna dimensions.', numel(A), Mtot_full);
end

g_full = reshape(g_tilde(1,1:S,:,1,d_idx), S, num_antennas);
G_R_full = permute(cat(2, real(g_full), imag(g_full)), [2,1]);   % Mtot_full x S

r_dim = Mtot_full - numel(A);
Q_A = null(G_R_full(:,A).').';   % r_dim x Mtot_full

breve_g_ret = Q_A * G_R_full(:,retained);   % r_dim x (S-|A|)

m_ret = mi_5d(1,retained,1,1,d_idx);
Dmat_ret = diag(m_ret.^2 * gamma_w);

B_A = breve_g_ret * Dmat_ret * breve_g_ret.';
B_A = (B_A + B_A.')/2;

[U_A, Lam_A] = eig(B_A/gamma_n);
[lam_sorted, idx] = sort(diag(Lam_A), 'descend');
U_A = U_A(:,idx);
lambda_vals_A = lam_sorted;
W_A = (1/sqrt(gamma_n)) * U_A.';

mu_ret = m_ret.^2;
WGmu_A = W_A * breve_g_ret * mu_ret.';   % r_dim x 1

Qn_cache = cell(r_dim,1);
for m_idx = 1:r_dim
    mag_sqr_H_m = Hm_arr(omega, lambda_vals_A(m_idx));
    [~, Qn_m] = get_time_domain(mag_sqr_H_m, dt, 1);
    Qn_cache{m_idx} = Qn_m;
end

geom.Q_A = Q_A;
geom.r_dim = r_dim;
geom.W_A = W_A;
geom.WGmu_A = WGmu_A;
geom.lambda_vals_A = lambda_vals_A;
geom.Qn_cache = Qn_cache;
geom.num_antennas = num_antennas;
end

function [alpha_hat, t0_hat, Lambda_hat] = estimate_with_geometry(geom, y, K, N, Tp, norm_fact, ...
    t, t0_true, mfTemplateFFT_raw, D_template, trial_idx, d_idx, dt)
% Applies a PRECOMPUTED geometry to ONE trial's data -- no eig, no null-space
% construction, no get_time_domain here. This is the only part redone per trial.

n_ret = numel(trial_idx);
r_dim = geom.r_dim;

y_dep = y(:,:,:,trial_idx,d_idx);
y_sq = reshape(y_dep, K, geom.num_antennas, n_ret);
y_R_full = permute(cat(2, real(y_sq), imag(y_sq)), [2 1 3]);   % Mtot_full x K x n_ret
y_R_A = pagemtimes(geom.Q_A, y_R_full);                          % r_dim x K x n_ret

z_A = pagemtimes(geom.W_A, y_R_A);   % r_dim x K x n_ret

mf_with_z_sum = zeros(1,1,K,n_ret);
for m_idx = 1:r_dim
    b_m = geom.WGmu_A(m_idx);
    Omega_t0 = reshape(b_m * D_template, K,1,K);
    z_m = reshape(z_A(m_idx,:,:), 1, K, 1, n_ret);
    mf_with_z_sum = mf_with_z_sum + dt*dt*pagemtimes(pagemtimes(z_m, reshape(geom.Qn_cache{m_idx},K,K,1,1)), Omega_t0);
end

[~, I] = max(mf_with_z_sum,[],3);
t0_hat = reshape((I-1)*dt, 1,1,n_ret);
t0_for_alpha = t0_true*ones(1,1,n_ret);

Af = fft(sensor_signal(t-t0_for_alpha, Tp, norm_fact), N, 1);
lag0 = round(Tp/dt);
Y_tensor = dt*ifft(Af .* mfTemplateFFT_raw, [], 1);
Rss_tensor = Y_tensor(lag0+1:(lag0+K),:,:);

num = zeros(1,1,n_ret);
denom = zeros(1,1,n_ret);
for m_idx = 1:r_dim
    Qn_m = geom.Qn_cache{m_idx};
    b_m = geom.WGmu_A(m_idx);
    Omega_m = b_m * Rss_tensor;
    resh_Omega = reshape(Omega_m,1,K,n_ret);
    z_m = reshape(z_A(m_idx,:,:), 1, K, n_ret);
    num = num + dt*dt*pagemtimes(pagemtimes(z_m,reshape(Qn_m,K,K,1)),pagetranspose(resh_Omega));
    denom = denom + dt*dt*pagemtimes(pagemtimes(resh_Omega,reshape(Qn_m,K,K,1)),pagetranspose(resh_Omega));
end
alpha_hat = num ./ denom;
Lambda_hat = (alpha_hat.^2) .* denom;
end

function [unique_sets, group_idx] = unique_cell_sets(set_cell)
% Groups a cell array of numeric row-vectors by exact content, regardless
% of length (MATLAB's built-in unique() doesn't handle this directly).
keys = cellfun(@(x) mat2str(sort(x)), set_cell, 'UniformOutput', false);
[unique_keys, ~, group_idx] = unique(keys);
unique_sets = cellfun(@(k) str2num(k), unique_keys, 'UniformOutput', false); %#ok<ST2NM>
end

function [alpha_hat, t0_hat, Lambda_hat] = null_and_estimate(A, y, mi_5d, g_tilde, gamma_w, gamma_n, ...
    Hm_arr, omega, dt, K, N, Tp, norm_fact, t, t0_true, mfTemplateFFT_raw, trial_idx, S, d_idx)
% Nulls the (possibly multi-sensor) set A -- no validation, applied at face value.

num_antennas = size(g_tilde,3);
Mtot_full = 2*num_antennas;
retained = setdiff(1:S, A);
n_ret = numel(trial_idx);

g_full = reshape(g_tilde(1,1:S,:,1,d_idx), S, num_antennas);
G_R_full = permute(cat(2, real(g_full), imag(g_full)), [2,1]);   % Mtot_full x S

% --- Q_A: null EVERY sensor in A simultaneously (generalizes directly to |A|>1) ---
r_dim = Mtot_full - numel(A);
Q_A = null(G_R_full(:,A).').';   % r_dim x Mtot_full, orthonormal rows

breve_g_ret = Q_A * G_R_full(:,retained);   % r_dim x (S-|A|)

y_dep = y(:,:,:,trial_idx,d_idx);
y_sq = reshape(y_dep, K, num_antennas, n_ret);
y_R_full = permute(cat(2, real(y_sq), imag(y_sq)), [2 1 3]);   % Mtot_full x K x n_ret
y_R_A = pagemtimes(Q_A, y_R_full);                              % r_dim x K x n_ret

m_ret = mi_5d(1,retained,1,1,d_idx);
Dmat_ret = diag(m_ret.^2 * gamma_w);

B_A = breve_g_ret * Dmat_ret * breve_g_ret.';
B_A = (B_A + B_A.')/2;

[U_A, Lam_A] = eig(B_A/gamma_n);
[lam_sorted, idx] = sort(diag(Lam_A), 'descend');
U_A = U_A(:,idx);
lambda_vals_A = lam_sorted;
W_A = (1/sqrt(gamma_n)) * U_A.';

z_A = pagemtimes(W_A, y_R_A);   % r_dim x K x n_ret

mu_ret = m_ret.^2;
WGmu_A = W_A * breve_g_ret * mu_ret.';   % r_dim x 1

t0_grid = reshape(t,1,1,[]);
[~, R00_full] = mf_integral_fft(sensor_signal(t-t0_grid,Tp,norm_fact), sensor_signal(t,Tp,norm_fact), 1, 1, K, dt, Tp);
R00_full = squeeze(R00_full);   % K x K

mf_with_z_sum = zeros(1,1,K,n_ret);
Qn_cache = cell(r_dim,1);
for m_idx = 1:r_dim
    mag_sqr_H_m = Hm_arr(omega, lambda_vals_A(m_idx));
    [~, Qn_m] = get_time_domain(mag_sqr_H_m, dt, 1);
    Qn_cache{m_idx} = Qn_m;
    b_m = WGmu_A(m_idx);
    Omega_t0 = reshape(b_m * R00_full, K,1,K);
    z_m = reshape(z_A(m_idx,:,:), 1, K, 1, n_ret);
    mf_with_z_sum = mf_with_z_sum + dt*dt*pagemtimes(pagemtimes(z_m, reshape(Qn_m,K,K,1,1)), Omega_t0);
end

[~, I] = max(mf_with_z_sum,[],3);
t0_hat = reshape((I-1)*dt, 1,1,n_ret);
t0_for_alpha = t0_true*ones(1,1,n_ret);

Af = fft(sensor_signal(t-t0_for_alpha, Tp, norm_fact), N, 1);
lag0 = round(Tp/dt);
Y_tensor = dt*ifft(Af .* mfTemplateFFT_raw, [], 1);
Rss_tensor = Y_tensor(lag0+1:(lag0+K),:,:);

num = zeros(1,1,n_ret);
denom = zeros(1,1,n_ret);
for m_idx = 1:r_dim
    Qn_m = Qn_cache{m_idx};
    b_m = WGmu_A(m_idx);
    Omega_m = b_m * Rss_tensor;
    resh_Omega = reshape(Omega_m,1,K,n_ret);
    z_m = reshape(z_A(m_idx,:,:), 1, K, n_ret);
    num = num + dt*dt*pagemtimes(pagemtimes(z_m,reshape(Qn_m,K,K,1)),pagetranspose(resh_Omega));
    denom = denom + dt*dt*pagemtimes(pagemtimes(resh_Omega,reshape(Qn_m,K,K,1)),pagetranspose(resh_Omega));
end
alpha_hat = num ./ denom;
Lambda_hat = (alpha_hat.^2) .* denom;
end