function [unique_sets, group_idx] = unique_cell_sets(set_cell)
% Groups a cell array of numeric row-vectors by exact content, regardless
% of length (MATLAB's built-in unique() doesn't handle this directly).
keys = cellfun(@(x) mat2str(sort(x)), set_cell, 'UniformOutput', false);
[unique_keys, ~, group_idx] = unique(keys);
unique_sets = cellfun(@(k) str2num(k), unique_keys, 'UniformOutput', false); %#ok<ST2NM>
end