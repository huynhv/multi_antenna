% GET_TIME_DOMAIN Returns the time-domain matrix of a given frequency domain
% vector.
function [Qn,Qn_matrix] = get_time_domain(mag_sqr_H, dt, nDeployments)
    mag_sqr_H_2 = [mag_sqr_H(1,1,:) mag_sqr_H(1,2:end,:) fliplr(conj(mag_sqr_H(1,2:end,:)))];
    Qn = ifft(mag_sqr_H_2/dt,[],2);
    L = floor(size(Qn,2)/2);
    mirrored_Qn = [Qn(1,L+2:end,:), Qn(1,1:L+1,:)];
    
    if nargout > 1
        toep_idx = (L+1) + (0:L) - (0:L).';   % Qn_matrix(i,j) = mirrored_Qn(L+1+j-i)
        Qn_matrix = reshape(mirrored_Qn(1, toep_idx(:), :), L+1, L+1, []);
    end
end