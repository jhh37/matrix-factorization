function [U, V, cost, iter] = xt_db_balm(M, W, r, U, V, varargin)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

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

[m, n] = size(M);

cost_best = inf;
error_best = inf;
M(W==0) = 0;
gamma = 5;
tau = 0.50;
rho = 50;
Lambda = zeros(m, r);

%% Initialization
Mz = U * V';
Mz(W~=0) = M(W~=0);

%% Optimization
% Outer loop
for iter = 1 : opts.max_iter
  
  % inner Gauss-Seidel iteration
  for trial = 1 : opts.max_trials
    % Step 1 : Update U_hat.
    U_hat = U - Lambda / rho;
    [U_hat, ~] = qr(U_hat, 0);
    
    % Step 2 : Update Q using Manopt.
    %     G = [ sqrt(rho / 2) * Mz' ; (U_hat + Lambda / rho)' ];
    %
    %     % Generate the problem data.
    %     GTG = G' * G;
    %
    %     % Create the problem structure.
    %     % manifold = spherefactory(m);
    %     manifold = stiefelfactory(m, r);
    %     problem.M = manifold;
    %
    %     % Define the problem cost function and its Euclidean gradient.
    %     problem.cost  = @(Q) - trace(Q' * (GTG * Q));
    %     problem.egrad = @(Q) -2 * GTG * Q;
    %
    %     % Numerically check gradient consistency (optional).
    %     % checkgradient(problem);
    %
    %     % Solve.
    %     tropts.verbosity = 0;
    %     [Q, Qcost, info, options] = trustregions(problem, [], tropts);
    
    % Step 2 : Update Q using SVD.
    G = [ sqrt(rho / 2) * Mz' ; (U_hat + Lambda / rho)' ];
    [~, ~, Q] = svd(G, 0);
    Q = Q(:, 1 : r);
    
    % Step 2 : Update U & V.
    S = G * Q;
    A = S(n+1 : end, :)';
    U = Q * A;
    V = sqrt(2 / rho) * S(1 : n, :) / A';
    
    % Step 3 : Update Z.
    Mz = U * V';
    Mz(W~=0) = M(W~=0);
  end
  
  error = norm(U - U_hat, 'fro') ^ 2;
  
  % PRIMAL-DUAL UPDATE
  if error <= tau * error_best
    % If constraint is pretty good, update the dual parameters
    Lambda = Lambda - rho * (U - U_hat);
    % rho = rho;
    error_best = error;
  else
    % Otherwise, increase the penalty proportion.
    rho = min(gamma * rho, 1e+20);
  end
  
  % STOPPING CRITERIA
  cost = sqrt(norm(W .* (M - U*V'), 'fro')^2 / sum(W(:)~=0));
  % disp(k); disp(cost);
  if opts.display, fprintf('[%03d] %0.6d\n', iter, cost);
  end
  if abs(cost - cost_best) < 1e-9, break
  else cost_best = cost;
  end
  
end

end