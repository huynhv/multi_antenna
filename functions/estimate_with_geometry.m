function [alpha_hat, t0_hat, Lambda_hat] = estimate_with_geometry(geom, y, K, N, Tp, norm_fact, ...
    t, t0_true, mfTemplateFFT_raw, D_template, trial_idx, d_idx, dt)
% Applies a PRECOMPUTED geometry to ONE trial's data -- no eig, no null-space
% construction, no get_time_domain here. This is the only part redone per trial.

n_ret = numel(trial_idx);
y_sq = reshape(y(:,:,:,trial_idx,d_idx), K, geom.num_antennas, n_ret);
y_R_full = permute(cat(2, real(y_sq), imag(y_sq)), [2 1 3]);   % Mtot_full x K x n_ret
z_A = pagemtimes(geom.W_A * geom.Q_A, y_R_full);               % r_dim x K x n_ret

% v = sum_m b_m*z_m*Qn_m, Q_sum = sum_m b_m^2*Qn_m
v = zeros(n_ret, K);
Q_sum = zeros(K, K);
for m_idx = 1:geom.r_dim
    b_m = geom.WGmu_A(m_idx);
    v = v + b_m * (reshape(z_A(m_idx,:,:), K, n_ret).' * geom.Qn_cache{m_idx});
    Q_sum = Q_sum + b_m^2 * geom.Qn_cache{m_idx};
end

[~, I] = max(v * D_template, [], 2);
t0_hat = reshape((I-1)*dt, 1, 1, n_ret);

% Alpha uses t0_true for every trial, so Rss is a single K x 1 vector.
Y = dt*ifft(fft(sensor_signal(t - t0_true, Tp, norm_fact), N, 1) .* mfTemplateFFT_raw);
lag0 = round(Tp/dt);
Rss = Y(lag0+1:lag0+K);
denom = dt*dt * (Rss.' * Q_sum * Rss);
alpha_hat  = reshape(dt*dt * (v * Rss) / denom, 1, 1, n_ret);
Lambda_hat = alpha_hat.^2 * denom;
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