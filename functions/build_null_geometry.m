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

% Same scaling as the main script: server noise is gamma_n/2 per real dimension.
[U_A, Lam_A] = eig(B_A);
[lam_sorted, idx] = sort(diag(Lam_A), 'descend');
lambda_vals_A = (2/gamma_n) * lam_sorted;
W_A = (1/sqrt(0.5*gamma_n)) * U_A(:,idx).';

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