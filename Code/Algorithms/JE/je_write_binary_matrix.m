function je_write_binary_matrix(filename, matrix, precision)
% JE_WRITE_BINARY_MATRIX
%   Writes a binary matrix in CSC binary format.

if nargin < 4, precision = 'double';
end

fid = fopen(filename, 'w');
fwrite(fid, matrix, precision);
fclose(fid);

end

