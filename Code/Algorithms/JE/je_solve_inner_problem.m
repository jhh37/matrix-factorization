function [X, AtQ, AtR] = je_solve_inner_problem(U, opts, substore)
% JE_SOLVE_LINEAR_LSQ
%   Solve the inner linear least squares problem.

% If the regularization parameter is 0, use the unregularized LLSQ solver.
if opts.nu == 0
  if opts.use_qr
    [X, AtQ, AtR] = je_solve_unreg_linear_lsq_qr(U, opts, substore);
  else
    [X, AtQ] = je_solve_unreg_linear_lsq(U, opts, substore);
  end
else
  [X, AtR] = je_solve_reg_linear_lsq(U, opts, substore);
end

end