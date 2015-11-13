function X = je_read_binary_matrix(filename, nRows, nCols, precision)
% JE_READ_BINARY_MATRIX
%   Read CSC binary matrix as matrix with specified precision. (default
%   precision: double)

if nargin < 4, precision = 'double';
end

fid = fopen(filename, 'r');
X = fread(fid, [nRows nCols], precision);
fclose(fid);

end

