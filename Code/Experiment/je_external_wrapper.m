function [U, V, cost, iter] = je_external_wrapper(M, W, r, U, varargin)
% JE_EXTERNAL_WRAPPER
%   An implementation of the Variable Projection algorithm on Low-rank
%   Matrix Completion Problem.

% tic
%% PRE-PROCESSING
% Fetch all the options.
opts = au_opts( ...
  'nu=0', ... % regularization parameter
  'max_iter=300', ... % maximum number of iterations
  'max_trials=50', ... % maximum number of trials
  'tol=1e-10', ... % function tolerance
  'display=1', ... % whether to display progress
  'alg=0', ... % which algorithm to use
  'r_init=0', ... % initial rank (for rank continuation-based approaches)
  varargin{:} );

% Fetch U if it has not been generated.
if size(U, 1) == 1 && size(U, 2) == 1
  if opts.r_init
    % If rank-continuation used, set U to be of rank r_init.
    U = je_generate_sample(U, size(M, 1), size(M, 2), opts.r_init);
  else
    % Otherwise, set U to be of rank r.
    U = je_generate_sample(U, size(M, 1), size(M, 2), r);
  end
end

%% ALGORITHM SELECTION
if strcmpi(opts.alg, 'PG_CSF')
  % Gotardo's CSF (Damped RW2): almost original
  M(W==0) = nan;
  [U, V, ~, cost, iter] = xt_pg_csf(M, r, U, opts.max_iter, opts.tol, opts.display);
elseif strcmpi(opts.alg, 'TO_DW')
  % Okatani's DW (Damped RW2 + constraint B): almost original
  % Tolerance is doubled as the script uses absolute cost value.
  [V, U, cost, iter] = xt_to_dw(M', W', r, U, 2 * opts.tol, opts.max_iter, opts.display);
  cost = sqrt(cost / sum(W(:) ~= 0));
elseif strcmpi(opts.alg, 'CH_LM_S')
  % Chen's LM_S (Damped RW0): modified
  % Tolerance is doubled as the script uses absolute cost value.
  [cost, iter, U] = xt_ch_lm_s(M', W', r, U, opts.max_iter, 2 * opts.tol, opts.display);
  cost = sqrt(cost / sum(W(:) ~= 0));
  V = NaN;
elseif strcmpi(opts.alg, 'CH_LM_S_GN')
  % Chen's LM_S_GN (Damped RW1): newly implemented based on CH_LM_S
  % Tolerance is doubled as the script uses absolute cost value.
  [cost, iter, U] = xt_ch_lm_s_gn(M', W', r, U, opts.max_iter, 2 * opts.tol, opts.display);
  cost = sqrt(cost / sum(W(:) ~= 0));
  V = NaN;
elseif strcmpi(opts.alg, 'CH_LM_S_RW2')
  % Chen's LM_S_RW2 (Damped RW2): newly implemented based on CH_LM_S
  % Tolerance is doubled as the script uses absolute cost value.
  [cost, iter, U] = xt_ch_lm_s_rw2(M', W', r, U, opts.max_iter, 2 * opts.tol, opts.display);
  cost = sqrt(cost(end) / sum(W(:) ~= 0));
  V = NaN;
elseif strcmpi(opts.alg, 'CH_LM_M')
  % Chen's LM_M (Damped reduced RW0): modified
  % Tolerance is doubled as the script uses absolute cost value.
  [cost, iter, U] = xt_ch_lm_m(M', W', r, U, opts.max_iter, 2 * opts.tol, opts.display);
  cost = sqrt(cost(end) / sum(W(:) ~= 0));
  V = NaN;
elseif strcmpi(opts.alg, 'CH_LM_M_GN')
  % Chen's LM_M_GN (Damped reduced RW1): modified
  % Tolerance is doubled as the script uses absolute cost value.
  [cost, iter, U] = xt_ch_lm_m_gn(M', W', r, U, opts.max_iter, 2 * opts.tol, opts.display);
  cost = sqrt(cost / sum(W(:) ~= 0));
  V = NaN;
elseif strcmpi(opts.alg, 'CH_LM_M_RW2')
  % Chen's LM_M_RW2 (Damped reduced RW2): newly implemented based on
  % CH_LM_M
  % Tolerance is doubled as the script uses absolute cost value.
  [cost, iter, U] = xt_ch_lm_m_rw2(M', W', r, U, opts.max_iter, 2 * opts.tol, opts.display);
  cost = sqrt(cost / sum(W(:) ~= 0));
  V = NaN;
elseif strcmpi(opts.alg, 'CO_LM_S')
  % Chen's LM_S (Damped RW0): original
  % CO_LM_S
  % Tolerance is doubled as the script uses absolute cost value.
  [U, ~] = qr(U, 0);
  [cost, iter, U] = xt_co_lm_s(M', W', r, U, opts.max_iter, 2 * opts.tol);
  cost = sqrt(cost(end) / sum(W(:) ~= 0));
  V = NaN;
elseif strcmpi(opts.alg, 'CO_LM_M')
  % Chen's LM_M (Damped reduced RW0): original
  % CO_LM_M
  % Tolerance is doubled as the script uses absolute cost value.
  [U, ~] = qr(U, 0);
  [cost, iter, U] = xt_co_lm_m(M', W', r, U, opts.max_iter, 2 * opts.tol);
  cost = sqrt(cost(end) / sum(W(:) ~= 0));
  V = NaN;
elseif strcmpi(opts.alg, 'CO_LM_M_GN')
  % Chen's LM_M_GN (Damped reduced RW1): original
  % CO_LM_M_GN
  % Tolerance is doubled as the script uses absolute cost value.
  [U, ~] = qr(U, 0);
  [cost, iter, U] = xt_co_lm_m_gn(M', W', r, U, opts.max_iter, 2 * opts.tol);
  cost = sqrt(cost(end) / sum(W(:) ~= 0));
  V = NaN;
elseif strcmpi(opts.alg, 'NB_RTRMC')
  % Boumal's RTRMC: cost function modified
  % NB_RTRMC
  % Tolerance is doubled as the script uses absolute cost value.
  [U, ~] = qr(U, 0);
  [U, V, cost, iter] = xt_nb_rtrmc_wrapper(M, W, r, U, ...
    ['max_iter=', num2str(opts.max_iter)], ...
    ['tol=', num2str(2 * opts.tol, '%e')], ...
    ['display=', num2str(opts.display)]);
  cost = sqrt(cost(end) / sum(W(:) ~= 0));
elseif strcmpi(opts.alg, 'RC_ALM')
  % Cabral's ALM: newly implemented based on author's description
  % RC_ALM
  % Tolerance is doubled as the script uses absolute cost value.
  M(W==0) = 100;
  [U, V] = je_alternation(M, W, r, U, nan, 'max_iter=0', 'display=0');
  [U, V, cost, iter] = xt_rc_alm(M, W, r, U, V, ...
    ['max_iter=', num2str(opts.max_iter)], ...
    ['max_trials=', num2str(opts.max_trials)], ...
    ['tol=', num2str(2 * opts.tol, '%e')], ...
    ['display=', num2str(opts.display)], ...
    ['nu=', num2str(opts.nu, '%e')]);
  cost = sqrt(cost / sum(W(:) ~= 0));
elseif strcmpi(opts.alg, 'RC_RCALM')
  % Cabral's ALM with rank continuation: newly implemented based on author's description
  % RC_RCALM
  % Tolerance is doubled as the script uses absolute cost value.
  [U, ~] = qr(U, 0);
  [U, V] = je_alternation(M, W, opts.r_init, U, nan, 'max_iter=0', 'display=0');
  
  [U, V, cost, iter] = xt_rc_rcalm(M, W, r, U, V, ...
    ['max_iter=', num2str(opts.max_iter)], ...
    ['max_trials=', num2str(opts.max_trials)], ...
    ['tol=', num2str(2 * opts.tol, '%e')], ...
    ['display=', num2str(opts.display)], ...
    ['nu=', num2str(opts.nu, '%e')]);
  cost = sqrt(cost / sum(W(:) ~= 0));
elseif strcmpi(opts.alg, 'DB_BALM')
  % Del Bue's BALM: newly implemented based on author's description
  % The projection matrix is replaced by QR decomposition.
  % DB_BALM
  % Tolerance is doubled as the script uses absolute cost value.
  [U, ~] = qr(U, 0);
  [U, V] = je_alternation(M, W, opts.r_init, U, nan, 'max_iter=0', 'display=0');
  [U, V, cost, iter] = xt_db_balm(M, W, r, U, V, ...
    ['max_iter=', num2str(opts.max_iter)], ...
    ['max_trials=', num2str(opts.max_trials)], ...
    ['tol=', num2str(2 * opts.tol, '%e')], ...
    ['display=', num2str(opts.display)], ...
    ['nu=', num2str(opts.nu, '%e')]);
elseif strcmpi(opts.alg, 'DO_BALM')
  % Del Bue's BALM: original
  % The projection matrix is replaced by QR decomposition.
  % DO_BALM
  [U, ~] = qr(U, 0);
  M(W==0) = 0;
  [U, V] = je_alternation(M, W, opts.r_init, U, nan, 'max_iter=0', 'display=0');
  [V, U, ~, ~, iter] = xt_do_balm_missing_bilinear_alter(M', W', V, U', 0, ...
    @xt_do_proj_qr, opts.max_iter, opts.max_trials);
  U = U';
  cost = sqrt(norm(W .* (M - U * V'), 'fro') ^ 2 / sum(W(:) ~= 0));
else
  error('No matching algorithm found.');
end

end

