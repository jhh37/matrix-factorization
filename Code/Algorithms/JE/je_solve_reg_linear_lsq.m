function [X, At, AtR] = je_solve_reg_linear_lsq(A, opts, substore)
% JE_SOLVE_REG_LINEAR_LSQ
%   Solve an unregularized linear least squares problem of the form: || Wt * (I kron A) * x - substore.mt ||^2_2

%% CHOLKESY DECOMPOSITION
r = size(A, 2);
At = cell(substore.unique.count, 1);
AtR = cell(substore.unique.count, 1);
if opts.use_indicator_weight
  for j = 1 : substore.unique.count
    At{j} = A(substore.unique.nz{j}, :);
    AtR{j} = At{j}' * At{j} + opts.nu * speye(r);
  end
else
  for j = 1 : substore.unique.count
    At{j} = bsxfun(@times, substore.unique.wnz{j}, A(substore.unique.nz{j}, :));
    AtR{j} = At{j}' * At{j} + opts.nu * speye(r);
  end
end

if opts.use_mp_inv
  AtP = cell(substore.unique.count, 1);
  for j = 1 : substore.unique.count
    AtP{j} = AtR{j} \ At{j}';
  end
end

%% NORMAL EQUATION
X = zeros(size(A, 2), substore.width);

if ~opts.use_mp_inv
  for i = 1 : substore.width
    j = substore.wid(i);
    X(:, i) = AtR{j} \ (At{j}' * substore.mt(substore.rng(i) : substore.rng(i + 1) - 1));
  end
else
  for i = 1 : substore.width
    j = substore.wid(i);
    X(:, i) = AtP{j} * substore.mt(substore.rng(i) : substore.rng(i + 1) - 1);
  end
end

X = X';

end