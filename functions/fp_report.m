function fp_report(fp_stats, verbosity_level)
% Print every entry in an FP-rate containers.Map, name-padded, via log_msg.
    log_msg(verbosity_level, 4, 'FP rate --');
    names = keys(fp_stats);
    for k = 1:numel(names)
        if names{k} ~= "combined"
            log_msg(verbosity_level, 4, '  %-14s %.4f', names{k}, fp_stats(names{k}));
        else
            continue
        end
    end
    log_msg(verbosity_level, 4, '  %-14s %.4f', "combined", fp_stats("combined"));
end