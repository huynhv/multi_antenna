function B = unsqueeze(A, dim)
%UNSQUEEZE Insert a singleton dimension into an array at the specified position.
%
%   B = UNSQUEEZE(A, DIM) inserts a singleton dimension (size 1) at
%   dimension DIM of array A.
%
%   Example:
%       A = rand(1,10,4);       % size: 1x10x4
%       B = unsqueeze(A, 3);    % size: 1x10x1x4

    sz = size(A);
    % pad size vector if dim is beyond current dimensions
    if dim > numel(sz) + 1
        sz(end+1:dim-1) = 1;
    end
    % insert singleton dimension
    newSize = [sz(1:dim-1) 1 sz(dim:end)];
    B = reshape(A, newSize);
end
