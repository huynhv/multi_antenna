% (REVISED MATCHED FILTERING METHOD) MF_INTEGRAL_FFT Returns the matched filter output of tempA and tempB using the FFT. 
function [sum_tensor, ui_tensor] = mf_integral_fft(A_tensor, B_tensor, ci, ai, K, dt, Tp)
    N = size(A_tensor,1) + size(B_tensor,1) - 1;
    Af = fft(A_tensor, N, 1);
    Bf = fft(B_tensor, N, 1);
    lag0 = round(Tp/dt);
    Y_tensor =  dt*ifft(Af .* Bf, [], 1);
    ui_tensor = Y_tensor(lag0+1 : (lag0+K), :, :, :, :); % need to make sure index is integer multiple
    sum_tensor = sum(ui_tensor .* ci .* ai, 2);
end