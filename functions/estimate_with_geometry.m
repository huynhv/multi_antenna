function [alpha_hat, t0_hat, Lambda_hat] = estimate_with_geometry(geom, y, K, N, Tp, norm_fact, ...
    t, t0_true, mfTemplateFFT_raw, D_template, trial_idx, d_idx, dt, offset_idx, alpha_uses_true_t0)
% Applies a PRECOMPUTED geometry to the selected trials -- no eig, no null-space
% construction, no get_time_domain here.
%   alpha_uses_true_t0 = false: joint estimation, alpha evaluated at each trial's t0_hat
%   alpha_uses_true_t0 = true : alpha evaluated at the true t0

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

% t0 search over the same window as the main estimator.
mf = v * D_template;                                 % n_ret x K
[~, I] = max(mf(:, offset_idx:end), [], 2);
col = I + offset_idx - 1;                            % n_ret x 1
t0_hat = reshape(t(col), 1, 1, n_ret);

if alpha_uses_true_t0
    Y = dt*ifft(fft(sensor_signal(t - t0_true, Tp, norm_fact), N, 1) .* mfTemplateFFT_raw);
    lag0 = round(Tp/dt);
    Rss = repmat(Y(lag0+1:lag0+K), 1, n_ret);        % K x n_ret, same for every trial
else
    Rss = D_template(:, col);                        % K x n_ret, rho(t - t0_hat) per trial
end
num   = dt*dt * sum(v.' .* Rss, 1);                  % 1 x n_ret
denom = dt*dt * sum(Rss .* (Q_sum * Rss), 1);        % 1 x n_ret
alpha_hat  = reshape(num ./ denom, 1, 1, n_ret);
Lambda_hat = alpha_hat.^2 .* reshape(denom, 1, 1, n_ret);
end