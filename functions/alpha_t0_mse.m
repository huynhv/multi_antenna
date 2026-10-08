function out = alpha_t0_mse(alpha_est, t0_est, alpha_true, t0_true)
% 2 x nDeployments: row 1 = alpha MSE, row 2 = t0 MSE (averaged over trials, dim 3).
out = [reshape(mean((alpha_est - alpha_true).^2, 3), 1, []);
       reshape(mean((t0_est    - t0_true).^2,    3), 1, [])];
end