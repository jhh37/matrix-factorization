function [U, V, cost, iter] = xt_rc_rcalm(M, W, r, U, V, varargin)
% XT_RC_ALM
%   Cabral's Augmented Lagrangian Method (ALM) with rank continuation
%   Solves || W .* (M - Z) ||_F ^ 2 + (rho / 2) * || Z - U * V' ||_F ^ 2
%   s.t. Z = U V' by relaxing the constraint.
%   Initially, rank(U) and rank(V) is set to be 3 higher than predicted.
%   Newly implemented based on the author's description.

%% Fetch option parameters.
opts = au_opts( ...
  'nu=0', ... % regularization parameter
  'max_iter=300', ... % maximum number of iterations
  'max_trials=50', ... % maximum number of trials for each iteration
  'tol=1e-10', ... % function tolerance
  'display=1', ... % whether to display progress
  'r_init=0', ... % initial rank
  varargin{:} );

%% Initialize parameters.
[m, n] = size(W);
M(W==0) = 100;

% Set initial rank.
if ~opts.r_init,
  r_init = min(r + 2, min(m, n));
end
[U, V, cost, iter] = xt_rc_alm(M, W, r_init, U, V, varargin{:});
Z = U * V';

%% Rank-reduce
for  rr = r_init - 1 : -1 : r
  [U, S, V] = svd(Z, 'econ');
  s = diag(S);
  S = diag(sparse(sqrt(s(1 : rr, 1))));
  U = U(:, 1 : rr) * S;
  V = V(:, 1 : rr) * S;
  % [U, ~] = qr(U, 0);
  [U, V, cost, iter_rr] = xt_rc_alm(M, W, rr, U, V, varargin{:});
  iter = iter + iter_rr;
  Z = U * V';
  
end

end

