% SENSOR_SIGNAL Computes the sensor signal defined in the paper.
function s = sensor_signal(t, Tp, norm_fact)
    B = pi/Tp;
    s = sin(B*t);
    s(t<0 | t>=Tp) = 0;
    s = s*norm_fact;
end
