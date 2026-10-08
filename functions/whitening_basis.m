function [UB_T, lambda_B] = whitening_basis(G_R, m, gamma_w)
% Eigenbasis of B = G_R*diag(gamma_w*m.^2)*G_R.' (sensor-noise covariance at the
% real-stacked antennas), eigenvalues sorted descending. G_R: Mtot x S x D, m: S x D.
[Mtot, S, D] = size(G_R);
B = pagemtimes(G_R .* reshape(gamma_w*m.^2, 1, S, D), pagetranspose(G_R));
B = (B + pagetranspose(B))/2;
[U, Lam] = pageeig(B);
lambda_B = zeros(Mtot, D);
UB_T = zeros(Mtot, Mtot, D);
for d = 1:D
    [lambda_B(:,d), idx] = sort(diag(Lam(:,:,d)), 'descend');
    UB_T(:,:,d) = U(:,idx,d).';
end
end