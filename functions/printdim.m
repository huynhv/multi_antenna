function printdim(x)
    name = inputname(1);
    fprintf('%s: %s\n', name, mat2str(size(x)));
end