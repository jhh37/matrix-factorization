function [U, V, cost, iter, store] = je_wiberg(M, W, r, U, store, varargin)
% JE_WIBERG
%   An implementation of the unregularized variable projection algorithm on low-rank
%   matrix completion problem

% tic
%% PRE-PROCESSING
% Custom function
vec = @(X) X(:);

% Fetch all the options.
opts = au_opts( ...
  'alg_type=2', ... % which VP algorithm to choose (1: Exact Wiberg, 2: Approximate Wiberg, 3: Approximate-approximate Wiberg)
  'retraction=1', ...
  'constraint=projection', ...
  'nu=0', ... % regularization parameter
  'lambda=1e-4', ... % damping parameter
  'lambda_inc_factor=10', ... % damping increase rate
  'lambda_dec_factor=10', ... % damping decrease rate
  'lambda_min=1e-14', ... % minimum value for lambda
  'max_iter=300', ... % maximum number of iterations
  'max_trials=50', ... % maximum number of trials for each iteration
  'tol=1e-10', ... % function tolerance
  'use_qr=1', ... % whether to use QR decomposition
  'use_mp_inv=0', ... % whether we are going to build explicit Moore-Penrose inverse matrices
  'use_indicator_weight=1', ... % whether we are going to do 0.5 ALS only
  'search_unique_weights=1', ... % whether to search for unique weight columns and rows.
  'use_layerwise_computation=1', ... % whether to compute Hessian by summing n layers or by using sparse matrix computation.
  'init_v_only=0', ... % whether we are going to do 0.5 ALS only
  'display=1', ... % whether to display progress
  varargin{:} );

% If store variable has not been passed, pre-process the required variables.
if nargin < 5 || ~isstruct(store)
  [opts, store] = je_preprocess_values(M, W, r, opts);
end

% Fetch U if it has not been generated.
if size(U, 1) == 1 && size(U, 2) == 1
  U = je_generate_sample(U, store.dim.m, store.dim.n, r);
end

% If manifold retraction is enabled, take the q-factor of U.
if opts.retraction, [U, ~] = qr(U, 0);
end
% toc

% Get the initial estimate of V.
% tic
[V, UtQ, UtR] = je_solve_inner_problem(U, opts, store.wc);
% toc
% Compute the cost: this will be replaced with a more efficient code later on.
[R, cost] = je_compute_error_matrix(M, W, U, V, store.dim.nnz);

%% SECOND-ORDER SUBPROBLEM SOLVER
% Initialize parameters.
iter = 0;
eval = 0;
lambda = opts.lambda;

% If display is on, show the initial cost.
if opts.display, fprintf('[%05d][%05d]\t%.2e\t%.6f\n', iter, eval, lambda, cost);
end

% Start iteration.
for iter = 1 : opts.max_iter
  
  % Compute the gradient and JTJ.
  if ~opts.use_indicator_weight, R = W .* R;
  end
  JTe = - vec((R * V)');
  JTJ = je_compute_unreg_wiberg_hessian(U, UtQ, UtR, V, R, r, opts, store);
  
  % Evaluate until the cost decreases or reaches the max number of fail.
  for trial = 1 : opts.max_trials
    eval = eval + 1;
    
    % Compute U*.
    dU = reshape((JTJ + lambda * speye(store.dim.mr)) \ JTe, r, store.dim.m)';
    U_eval = U + dU;
    
    % If manifold retraction is enabled, take the q-factor of U + dU.
    if opts.retraction, [U_eval, ~] = qr(U_eval, 0);
    end
    
    % U_eval = reshape((JTJ + lambda * speye(store.dim.mr)) \ (JTe + lambda * vec(U')), r, store.dim.m)';
    
    % Compute the corresponding V*.
    [V_eval, UtQ, UtR] = je_solve_inner_problem(U_eval, opts, store.wc);
    
    % Compute the cost: this will be replaced with a more efficient code later on.
    [R, cost_eval] = je_compute_error_matrix(M, W, U_eval, V_eval, store.dim.nnz);
    
    if cost_eval < cost, break
    else
      lambda = lambda * opts.lambda_inc_factor;
    end
  end
  
  % If display is on, print the number of iterations and evaluations and the cost function value.
  cost_diff = abs(cost - cost_eval);
  if opts.display, fprintf('[%05d][%05d] %.2e %.6f %.6e\n', iter, eval, lambda, cost_eval, cost_diff);
  end
  
  % Update variables.
  U = U_eval;
  V = V_eval;
  lambda = max(lambda / opts.lambda_dec_factor, opts.lambda_min);
  less_than_tol = cost_diff < cost * opts.tol;
  cost = cost_eval;
  
  % Stop if the change in cost function is smaller than the function
  % tolerance value.
  if trial == opts.max_trials || less_than_tol, break;
  end
end

end

