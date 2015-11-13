function [U, V, cost, iter, store] = je_alternation(M, W, r, U, store, varargin)
% JE_ALTERNATION
%   An implementation of the alternation algorithm optimized for medium-sized matrices.
%   Mainly exploits the structure of the matrices and repetitions in arithmetics.

%% PRE-PROCESSING

% Fetch all the options.
opts = au_opts( ...
  'retraction=1', ... % q factor-based manifold retraction
  'nu=0', ... % regularization parameter
  'max_iter=300', ... % maximum iteration
  'tol=1e-10', ... % function tolerance
  'use_qr=1', ... % whether to use QR decomposition
  'use_mp_inv=0', ... % whether we are going to build explicit Moore-Penrose inverse matrices
  'use_indicator_weight=1', ... % whether we are going to do 0.5 ALS only
  'use_layerwise_computation=0', ... % whether to compute Hessian by summing n layers or by using sparse matrix computation.
  'search_unique_weights=1', ... % whether to search for unique weight columns and rows.
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
V = je_solve_inner_problem(U, opts, store.wc);
% Compute the cost: this will be replaced with a more efficient code later on.
[~, cost] = je_compute_error_matrix(M, W, U, V, store.dim.nnz);
iter = 0;

% If display is on, show the cost.
if opts.display, fprintf('[%05d]\t%.6f\n', iter, cost);
end

if opts.init_v_only, return
end

%% ALTERNATION
for iter = 1:opts.max_iter
  
  U = je_solve_inner_problem(V, opts, store.wr);
  % If manifold retraction is enabled, take the q-factor of U.
  if opts.retraction, [U, ~] = qr(U, 0);
  end
  V = je_solve_inner_problem(U, opts, store.wc);
  
  % Compute the cost: this will be replaced with a more efficient code later on.
  [~, cost_eval] = je_compute_error_matrix(M, W, U, V, store.dim.nnz);
  % cost_eval = sqrt(norm(je_compute_errfunc(U, V, store.wc, opts)) ^ 2 / store.dim.nnz);
  
  % If display is on, show the cost.
  if opts.display, fprintf('[%05d]\t%.6f\n', iter, cost_eval);
  end
  
  % Stop if the convergence criteria is satisfied.
  if abs(cost - cost_eval) < cost * opts.tol
    break
  else cost = cost_eval;
  end
end

end

