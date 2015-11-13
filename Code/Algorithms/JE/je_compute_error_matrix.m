function [R, cost] = je_compute_error_matrix(M, W, A, B, nnz)
%JE_COMPUTE_ERRFUNC
%   Compute the error function such that:
%   ||f||^2 = || W .* (A * B' - M) ||^2_F

if nargin < 5, nnz =  sum(W(:) ~= 0);
end

%% ERROR MATRIX
R = W .* (A * B' - M);

if nargout > 1, cost = sqrt(norm(R, 'fro')^2 / nnz);
end

end
