function geom = cached_null_geometry(cache, A, mi_5d, g_tilde, gamma_w, gamma_n, Hm_arr, omega, dt, K, S, d_idx)
% Returns build_null_geometry(A, ...) for deployment d_idx, building it only on first request.
% CACHE is a containers.Map (a handle object), so entries added here persist for the caller.
% Only valid while gamma_w/gamma_n are fixed -- create a fresh cache per channel-SNR point.
key = sprintf('%d_%s', d_idx, mat2str(sort(A)));
if ~isKey(cache, key)
    cache(key) = build_null_geometry(A, mi_5d, g_tilde, gamma_w, gamma_n, Hm_arr, omega, dt, K, S, d_idx);
end
geom = cache(key);
end