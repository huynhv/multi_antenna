% Clear all variables and close all existing figures.
clearvars
close all

addpath('.\functions')

% Set seed for reproducibility
seed = 22;
rng(seed);

nDeployments = 10;
% Set the number of measurements per sensor, starting from 0
K = 251;
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
alpha_true = 2;
t0_true = 1;

% Define range for distribution of mi and ti (tau_i)
max_mi = 1;
min_mi = 0.5;

max_ti = 1;
min_ti = 0;

% Set maximum value of t0 such that the signal from each sensor is guaranteed to
% be in the interior of the observation period.
t0_max = T0-Tp-max_ti;

if t0_true > t0_max
    error("True value of t0 exceeds maximum!")
end

% Set the number of antennas
num_antennas = 8;
Mtot = 2*num_antennas;

%%% Always generate w, n, and channel gain for maximum S, then use subsets for different S trials
S_max = max(sensor_vals);
% Generate sensor noise with dimensions: K x S x trials x deployments.
all_w = randn(K,S_max,1,nTrials,nDeployments);
% Generate server noise.
n = (randn(K,1,num_antennas,nTrials,nDeployments) + 1i*randn(K,1,num_antennas,nTrials,nDeployments));
% Generate channel gains.
all_gi = (randn(1,S_max,num_antennas,1,nDeployments) + 1i*randn(1,S_max,num_antennas,1,nDeployments))/sqrt(2);

% Generate mi and ti.
all_mi = sort(unifrnd(min_mi,max_mi,1,S_max,1,1,nDeployments),2,'descend');
all_ti = sort(unifrnd(min_ti,max_ti,1,S_max,1,1,nDeployments),2,'ascend');

selected_schemes = "EPC";

% Define Rayleigh distribution parameters.
rayleigh_factor = 1/sqrt(2);
E_mag_g_sqr = 2*rayleigh_factor^2;

% Define expected value of mi^2.
E_mi_sqr = (1/12) * (max_mi-min_mi)^2 + 0.5*(min_mi+max_mi);

% Set miscellaneous parameters
norm_fact = 1;
epsilon = 1e-4;
Bee = pi/Tp;
w_psd_constant = 1;
n_psd_constant = 1;
N = 2*K - 1;

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

alpha_branch_fisher = @(bm, lambda) integral(@(w) bm.^2 .* mag_sqr_S0(w).^2 ./ Vm_arr(w, lambda), -Inf, Inf, 'ArrayValued', true);
t0_branch_fisher = @(bm, lambda) integral(@(w) bm.^2 .* (w .* mag_sqr_S0(w)).^2 ./ Vm_arr(w, lambda), -Inf, Inf, 'ArrayValued', true);

%% Define experiment list
experiment_list = "cwe_disjoint_peak_power_opt";

for experiment_idx = 1:numel(experiment_list)
    experiment = experiment_list(experiment_idx);

    pad_str = '----------';
    msg = [pad_str  ' ' 'Running: ' char(experiment) ' ' pad_str];
    fprintf('\n%s\n', msg);

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
    rho_empirical_var = rho_crlb;
    rho_empirical_mse = rho_crlb;
    rho_empirical_bias = rho_crlb;

    % Start run timer
    disp("Starting runtime...")
    loopTic = tic;
    
    % Iterate through agent snr values.
    for agent_db_idx = 1:length(agent_db_values)
        agent_db = agent_db_values(agent_db_idx);
    
        % Set sensor noise power. Note the factor of 1/dt, which is equivalent to
        % applying an anti-aliasing filter.
        gamma_w = (dt/w_psd_constant) * (Ps / db2magTen(agent_db));
        scaled_w = sqrt(gamma_w * w_psd_constant /dt) * all_w; % --> variance of this should be gamma_w/dt
        scaled_w_psd_constant = gamma_w * w_psd_constant; %% = Nw/2
        
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

                    %% Compute values for multi-antenna
                    m = reshape(mi_5d, S, nDeployments);           % S x nDeployments
                    g = reshape(g_tilde, S, num_antennas, nDeployments);       % S x num_antennas x nDeployments (complex, per-antenna compensated gain)

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

                    % Iterate through dropout values.
                    for dropout_idx = 1:length(dropout_vals)
                        save_dim = sensor_idx;

                        % Iterate through channel snr values.
                        for channel_db_idx = 1:size(channel_db_values,2)
                            channel_db_start = tic;
    
                            % Print statement for at-a-glance performance.
                            disp('')
                            disp("=== " + scheme + " ===")
                            disp("** Agent SNR = " + (agent_db_values(agent_db_idx)) + " dB **")
                            disp("** Channel SNR = " + (channel_db_values(scheme_idx,channel_db_idx,agent_db_idx)) + " dB **")
                            disp("** S = " + S + ", Dropout = " + dropout_vals(dropout_idx) + " **")

                            if scheme == "EPC"
                                gamma_n = (dt/n_psd_constant) * Ps * E_mi_sqr * E_mag_g_sqr / db2magTen(channel_db_values(scheme_idx,channel_db_idx,agent_db_idx));
                            end

                            % Set server noise power. Note the factor of 1/dt which is
                            % equivalent to applying an anti-aliasing filter.
                            scaled_n = sqrt( (gamma_n*n_psd_constant/dt) / 2 ) * n;
                            scaled_n_psd_constant = gamma_n * n_psd_constant; % == N0/2

                            %%%%%%%%%%%%%%%%% ESTIMATION %%%%%%%%%%%%%%%%%
                            total_est_start = tic;
    
                            disp("Begin ML estimation...")
                            for strat_idx = 1:num_strats
                                strat = all_strats(strat_idx);
                                disp("Running ML estimation for " + strat + " | obj = " + obj_func)
        
                                ai = 1;
                                bi = 0;

                                y = pagetranspose(out_cws + scaled_n);

                                %%% Baseline single-antenna case from
                                %%% previous manuscripts
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

                                    y_sq = reshape(y, K, num_antennas, nTrials, nDeployments);   % K x num_antennas x nTrials x nDeployments
                                    g_sq = reshape(g_tilde, S, num_antennas, nDeployments);   % S x num_antennas x nDeployments

                                    y_R = cat(2, real(y_sq), imag(y_sq));   % K x 2M x nTrials x nDeployments
                                    y_R = permute(y_R, [2 1 3 4]);   % 2M x K x nTrials x nDeployments
                                    
                                    G_R = cat(2, real(g_sq), imag(g_sq));
                                    G_R = permute(G_R, [2,1,3]);

                                    W_bcast = reshape(W, Mtot, Mtot, 1, nDeployments);  % broadcast over nTrials

                                    z = pagemtimes(W_bcast, y_R);   % 2M x K x nTrials x nDeployments, decorrelated signal+noise

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

                                    % Get time domain matrix for Qn.
                                    [max_y_vals,I] = max(mf_with_z_sum, [], 3);

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
                                end

                                % Compute empirical variance.
                                rho_empirical_var(strat_idx,agent_db_idx,channel_db_idx,1,scheme_idx,pivot_idx,save_dim,:) = var(alpha_estimates,0,3);
                                rho_empirical_var(strat_idx,agent_db_idx,channel_db_idx,2,scheme_idx,pivot_idx,save_dim,:) = var(t0_estimates_for_plot,0,3);

                                % Compute empirical mse.
                                rho_empirical_mse(strat_idx,agent_db_idx,channel_db_idx,1,scheme_idx,pivot_idx,save_dim,:) = mean((alpha_estimates - alpha_true).^2,3);
                                rho_empirical_mse(strat_idx,agent_db_idx,channel_db_idx,2,scheme_idx,pivot_idx,save_dim,:) = mean((t0_estimates_for_plot - t0_true).^2,3);

                                % Compute empirical bias.
                                rho_empirical_bias(strat_idx,agent_db_idx,channel_db_idx,1,scheme_idx,pivot_idx,save_dim,:) = mean(alpha_estimates,3);
                                rho_empirical_bias(strat_idx,agent_db_idx,channel_db_idx,2,scheme_idx,pivot_idx,save_dim,:) = mean(t0_estimates_for_plot,3);

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
                                end
                            end
                            total_est_time = toc(total_est_start);
                            disp("Total ML estimation time: " + floor(total_est_time/60) + " minutes " + mod(total_est_time,60) + " seconds")
                        end
                        total_channel_db_time = toc(channel_db_start);
                        disp("Total runtime for channel SNR = " + (channel_db_values(scheme_idx,channel_db_idx,agent_db_idx)) + " dB: " + floor(total_channel_db_time/60) + " minutes " + mod(total_channel_db_time,60) + " seconds")
                        fprintf("\n");
                    end
                end
            end
        end
    end
    
    % Stop run timer.
    runtime = toc(loopTic);
    % Display total run time.
    disp(experiment + " took " + floor(runtime/60) + " minutes " + mod(runtime,60) + " seconds")
    
    %% Compute averages across deployments.
    if nDeployments == 1
        avg_dep_var = rho_empirical_var;
        avg_dep_mse = rho_empirical_mse;
        avg_dep_crlb = rho_crlb;
        avg_dep_bias = rho_empirical_bias;
    else
        avg_dep_var = mean(rho_empirical_var,length(size(rho_empirical_var)));
        avg_dep_mse = mean(rho_empirical_mse,length(size(rho_empirical_mse)));
        avg_dep_crlb = mean(rho_crlb,length(size(rho_crlb)));
        avg_dep_bias = mean(rho_empirical_bias,length(size(rho_empirical_bias)));
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

            % Collect all plot values for y-axis scaling.
            all_vals = [avg_dep_crlb(selected_strat_idxs,agent_db_idx,:,param_idx,:,1:sensor_dimension);
                        avg_dep_var(selected_strat_idxs,agent_db_idx,:,param_idx,:,1:sensor_dimension);
                        avg_dep_mse(selected_strat_idxs,agent_db_idx,:,param_idx,:,1:sensor_dimension)
                        ];

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
                            plot(ax, x_axis_series,(squeeze(avg_dep_var(strat_idx,agent_db_idx,:,param_idx,scheme_idx,count_idx,1))),'-^','color',colorMap(iter_key),'LineWidth', plot_line_width)
                            plot(ax, x_axis_series,(squeeze(avg_dep_crlb(strat_idx,agent_db_idx,:,param_idx,scheme_idx,count_idx,1))),'--o','color',colorMap(iter_key),'LineWidth', plot_line_width)
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