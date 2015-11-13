function [U, V] = je_generate_sample(sample_id, m, n, r, base)
% JE_GENERATE_SAMPLE
%   Generate a set of samples with input dimension.
%

if nargin < 5, base = 1348032014;
end

width = 10;

rng(base + sample_id * width, 'twister');

U = randn(m, r);

if nargout > 1, V = randn(n, r);
end

end

