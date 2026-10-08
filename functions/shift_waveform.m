function y = shift_waveform(x, t, delta_t)
% SHIFT_WAVEFORM  Time-shift a waveform and zero out any part that spills
%                 outside the original observation window [t(1), t(end)].
%
%   y = shift_waveform(x, t, delta_t)
%
%   x       : samples of x(t), row or column vector
%   t       : time grid x was sampled on (assumed [0, T0], any spacing)
%   delta_t : shift amount. Positive delta_t DELAYS the pulse (moves it
%             later in time); negative ADVANCES it (moves it earlier).
%
%   y(t) = x(t - delta_t), but any t - delta_t that falls outside
%   [t(1), t(end)] (i.e. would require t<0 or t>T0 samples of x that
%   don't exist) is set to zero rather than extrapolated.

    t_query = t - delta_t;                      % where to sample x to get x(t - delta_t)
    y = interp1(t, x, t_query, 'linear', 0);     % 0 = fill value outside range -> zero

    % preserve original shape (row/col)
    y = reshape(y, size(x));
end