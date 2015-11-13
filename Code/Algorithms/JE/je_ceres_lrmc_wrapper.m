function [U, V, cost, iter, store] = je_ceres_lrmc_wrapper(M, W, r, U, store, varargin)
% JE_CERES_WRAPPER
%   The wrapper function for CERES-solver LMs

% tic
%% PRE-PROCESSING
% Fetch all the options.
opts = au_opts( ...
  'nu=0.0', ... % regularization parameter
  'retraction=0', ... % manifold retraction
  'max_iter=300', ... % maximum number of iterations
  'tol=1e-10', ... % function tolerance
  'use_qr=1', ... % whether to use QR decomposition
  'use_mp_inv=0', ... % whether we are going to build explicit Moore-Penrose inverse matrices
  'use_indicator_weight=1', ... % whether we are going to do 0.5 ALS only
  'use_layerwise_computation=0', ... % whether to compute Hessian by summing n layers or by using sparse matrix computation
  'search_unique_weights=1', ... % whether to search for unique weight columns and rows
  'init_v_only=1', ... % whether we are going to do 0.5 ALS only
  'max_als=0', ... % set the maximum number of additional set of ALS iterations
  'use_inner_iters=0', ... % whether to use Ceres inner iterations
  'use_auto_diff=0', ... % whether to use Ceres automatic differentiation
  'use_pca=0', ... % whether to use PCA formulation
  'reg_unreg=0', ... % whether to solve the regularized problem then the unregularized problem
  'sample_id=0', ... % whether to use the already-generated sample file
  'init_rank=0', ... % initial rank to start with (0 means start from the specified rank)
  'write_binary_mw=1', ... % whether to write binary files for the measurement matrix and the weight matrix
  'start_from_binary_uv=1', ... % whether to start directly from binary files of U and V
  'display=1', ... % whether to display progress
  varargin{:} );

% If the maximum number of ALS iterations is greater than 0, set
% init_v_only flag to 0.
if opts.max_als > 0, opts.init_v_only = 0;
end

% If store variable has not been passed, pre-process the required variables.
if nargin < 5 || ~isstruct(store) || opts.init_v_only
  [opts, store] = je_preprocess_values(M, W, r, opts);
end

% Fetch U if it has not been generated.
if size(U, 1) == 1 && size(U, 2) == 1
  U = je_generate_sample(U, store.dim.m, store.dim.n, r + opts.use_pca);
end

% If manifold retraction is enabled, take the q-factor of U.
if opts.retraction, [U, ~] = qr(U, 0);
end
% toc

% Get the initial estimate of V or ALS output of U and V.
% tic
if opts.init_v_only,
  [V, ~, ~] = je_solve_inner_problem(U, opts, store.wc);
  % Compute the cost: this will be replaced with a more efficient code later on.
  [~, cost] = je_compute_error_matrix(M, W, U, V, store.dim.nnz);
  iter = 0;
else
% tic
  [U, V, cost, iter] = je_alternation(M, W, r, U, nan, varargin{:}, ...
    sprintf('max_iter=%d', opts.max_als), ...
    sprintf('retraction=%d', opts.retraction));
  % If display is on, show the cost.
  if opts.display, fprintf('--- [ %d ALS iterations completed ] ---\n', iter);
  end
end
% toc

% If display is on, show the cost.
if opts.display, fprintf('[%05d]\t%.6f\n', 0, cost);
end

% If the sample ID is not given, write the sample matrices in CSC-binary format.
if opts.sample_id == 0, opts.sample_id = sprintf('s%08d', floor(rand * 1e8));
end

filename = ['Temp/', opts.sample_id, '_r', num2str(r)];

% Write the sample matrices in CSC-binary format.
if opts.write_binary_mw
  je_write_binary_matrix([filename, '_M.bin'], M);
  je_write_binary_matrix([filename, '_W.bin'], W);
end
je_write_binary_matrix([filename, '_U0.bin'], U);
je_write_binary_matrix([filename, '_V0.bin'], V);

% If rank-1 initialization is set, then run first with rank-1 constraint.
if opts.init_rank ~= 0
  U_init = je_ceres_wrapper(M, W, opts.init_rank, U(:, 1:(opts.init_rank + opts.use_pca)), tol, max_iter, varargin{:}, 'init_rank=0');
  U(:, 1 : opts.init_rank) = U_init(:, 1 : opts.init_rank);
  
  % If using PCA, also copy the translational vector.
  if opts.use_pca, U(:, end) = U_init(:, end);
  end
end

statement = ['./je_ceres_lrmc ',...
  '--dataset=', opts.sample_id, ' ',...
  '-m=', num2str(store.dim.m), ' ',...
  '-n=', num2str(store.dim.n), ' ',...
  '-r=', num2str(r), ' ',...
  '--tol=', num2str(opts.tol), ' ',...
  '--max_iter=', num2str(opts.max_iter), ' ',...
  '--nu=', num2str(opts.nu), ' ',...
  '--use-pca=', num2str(opts.use_pca), ' ',...
  '--use-inner-iters=', num2str(opts.use_inner_iters'), ' ',...
  '--use-auto-diff=', num2str(opts.use_auto_diff), ' ',...
  '--display=', num2str(opts.display), ' ',...
  '--debug=' num2str(opts.display), ' ',...
  ];

status = system(statement);

% Throw an error if not run correctly.
if status ~= 0, error('JE_CERES failed running');
end

% Read U and V.
U = je_read_binary_matrix([filename, '_U.bin'], store.dim.m, r + opts.use_pca);

% Increment the iteration counter.
iter = iter + je_read_binary_matrix([filename, '_iters.bin'], 2, 1);

% If using RU-mode, set current U and V to be the initial matrices for the
% un-regularized problem (warm start).
if opts.reg_unreg && (opts.nu > 0.0)
  [U, V, cost, iter_unreg] = je_ceres_lrmc_wrapper(M, W, r, U, store, varargin{:}, ...
    'reg_unreg=0', 'max_als=0', 'nu=0.0', ['sample_id=', opts.sample_id], 'write_binary_mw=0');
  iter = iter + iter_unreg;
else
  % Otherwise, fetch relevant statistics.
  V = je_read_binary_matrix([filename, '_V.bin'], store.dim.n, r);
  
  % If using PCA formulation, concatenate 1-vector to V.
  if opts.use_pca, V = [V ones(store.dim.n, 1)];
  end
  
  % Compute the final cost.
  [~, cost] = je_compute_error_matrix(M, W, U, V, store.dim.nnz);
  if opts.display, fprintf('[%05d]\t%.6f\n', iter, cost);
  end
  
  delete([filename, '_U.bin']);
  delete([filename, '_V.bin']);
  delete([filename, '_U0.bin']);
  delete([filename, '_V0.bin']);
  
  if opts.write_binary_mw
    delete([filename, '_M.bin']);
    delete([filename, '_W.bin']);
  end
  
end