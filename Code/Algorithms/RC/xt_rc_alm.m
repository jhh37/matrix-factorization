function [U, V, cost, iter] = xt_rc_alm(M, W, r, U, V, varargin)
% XT_RC_ALM
%   Cabral's Augmented Lagrangian Method (ALM)
%   Solves || W .* (M - Z) ||_F ^ 2 + (rho / 2) * || Z - U * V' ||_F ^ 2
%   s.t. Z = U V' by relaxing the constraint.
%   Newly implemented based on the author's description.

%% Fetch option parameters.
opts = au_opts( ...
  'nu=0', ... % regularization parameter
  'max_iter=300', ... % maximum number of iterations
  'max_trials=50', ... % maximum number of trials for each iteration
  'tol=1e-10', ... % function tolerance
  'display=1', ... % whether to display progress
  varargin{:} );

%% Initialize parameters.
[m, n] = size(M);
cost_best = inf;
error_best = inf;
widx = W(:)~=0;
M(~widx) = 100;
gamma = 1.01; % penalty increasing proportion
tau = 0.5; % minimum required error drop
rho = 0.5;  % initial penalty parameter
Y = zeros(m, n); % initial Lagrangian dual variables
Z = U * V'; % initially the constraint is satisfied.

%% Algorithm body
for iter = 1 : opts.max_iter
  
  % Update U, V and Z until || Z - U V' ||_F ^ 2 reaches some criteria.
  for trial = 1 : opts.max_trials
    
    Z2 = Z + Y / rho;
    if opts.nu ~= 0,
      U = Z2 * V / (V' * V + (opts.nu / rho) * speye(r));
      V = Z2' * U / (U' * U + (opts.nu / rho) * speye(r));
    else
      U = Z2 / V';
      V = (U \ Z2)';
    end
    X = U * V' - Y / rho;
    
    % Z(~widx) = X(~widx);
    % Z(widx) = (2 * M(widx) + rho * X(widx)) / (2 + rho);
    Z = W .* ((2 + rho) \ (2 * M + rho * X)) + (~W) .* X;
    
    % Stopping criteria: check the size of || Z - U V' ||_F ^ 2
    R = Z - U * V';
    error = norm(R(:), 'fro') ^ 2;
    if error <= tau * error_best
      error_best = error;
      break
    end
  end
  
  % Update dual variables.
  Y = Y + rho * R;
  % Lambda = max(zeros(m, n), min(Lambda + rho * R, 1e+20 * ones(m, n)));
  
  % Increase the penalty proportion.
  rho = min(gamma * rho, 1e+20);
  
  % Compute cost.
  cost = norm(W .* (M - U * V'), 'fro') ^ 2;
  if opts.display,
    fprintf('[%04d] %.6d (rank: %02d)\n', iter, sqrt(cost / sum(W(:) ~= 0)), r);
  end
  
  % Stopping criteria: check the cost
  cost_alm = cost + trace(Y' * (Z - U * V')) + rho / 2 * error;
  if opts.nu,
    cost_alm = cost_alm + 0.5 * opts.nu * (norm(U, 'fro') ^ 2 + norm(V, 'fro') ^ 2);
  end
  if abs(cost_alm - cost_best) < opts.tol, break
  else cost_best = cost_alm;
  end
end

end