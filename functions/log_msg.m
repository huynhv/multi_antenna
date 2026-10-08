function log_msg(verbosity_level, level, fmt, varargin)
% level 1 = major section (scheme/agent-SNR/experiment boundaries)
% level 2 = sub-step (channel SNR point, per-block timing)
% level 3 = detail/diagnostic (per-strategy, T_i/geometry/detection stats)
% level 4 = fine-grained (per-sensor, per-candidate detail)
% level 5 = trace (innermost loop, per-trial/per-iteration detail)
if level > verbosity_level
    return
end
indent = repmat('  ', 1, level-1);
switch level
    case 1
        prefix = sprintf('\n%s=== ', indent);
        suffix = ' ===';
    case 2
        prefix = sprintf('%s-- ', indent);
        suffix = '';
    case 3
        prefix = sprintf('%s.. ', indent);
        suffix = '';
    case 4
        prefix = sprintf('%s.... ', indent);
        suffix = '';
    case 5
        prefix = sprintf('%s...... ', indent);
        suffix = '';
    otherwise
        prefix = sprintf('%s   ', indent);
        suffix = '';
end
fprintf([prefix fmt suffix '\n'], varargin{:});
end