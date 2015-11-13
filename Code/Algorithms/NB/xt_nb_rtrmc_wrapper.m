function [U, V, cost, iter] = xt_nb_rtrmc_wrapper(M, W, r, U, varargin)

opts = au_opts( ...
  'nu=0', ... % regularization parameter
  'max_iter=300', ... % maximum number of iterations
  'max_trials=50', ... % maximum number of trials for each iteration
  'tol=1e-10', ... % function tolerance
  'display=1', ... % whether to display progress
  varargin{:} );


% Test code for the RTRMC/RCGMC algorithm (low-rank matrix completion).
% Algorithm by Nicolas Boumal and P.-A. Absil. Code by Nicolas Boumal,
% UCLouvain, April 2, 2013.
%
% http://perso.uclouvain.be/nicolas.boumal/RTRMC/
% http://sites.uclouvain.be/absil/RTRMC/
%
% See also: rtrmc buildproblem initialguess

% clear all;
% close all;
% clc;

% if exist('manopt_version', 'file') ~= 2
%   cd '../../';
%   
%   importmanopt;
%   cd ..;
%   warning('manopt:version', ['Using Manopt 1.0.7. This is probably not ' ...
%     'the most up-to-date version of Manopt. Please consider going to ' ...
%     'http://www.manopt.org to obtain the latest version.\n\n' ...
%     'Manopt 1.0.7 was added to your Matlab path.\n\n']);
% end

% If this script fails, try executing installrtrmc.m, to compile the mex
% files. This will not be necessary if launching the present script works
% out of the box, which there is a decent chance it should.

%% Problem instance generation

% Dimensions of the test problem
% load ~/Dropbox/Research/Projects/Beyond' Wiberg'/Datasets/Dino_trimmed/dino_trimmed.mat; r = 4;
% load ~/Dropbox/Research/Projects/Beyond'
% Wiberg'/Datasets/Giraffe/giraffe.mat; r = 6;

[m, n] = size(M);
% k = sum(W(:));

% Fetch U if it has not been generated.
if size(U, 1) == 1 && size(U, 2) == 1
  [U, ~] = qr(je_generate_sample(U, size(M, 1), size(M, 2), r), 0);
end

% m = 5000;                               % number of rows n = 5000;
% % number of columns r = 10;                                 % rank k =
% 4*r*(m+n-r);                        % number of known entries

% Generate an m-by-n matrix of rank true_rank in factored form: A*B
% true_rank = r;
% A = randn(m, true_rank)/true_rank.^.25; B = randn(true_rank,
% n)/true_rank.^.25;

% Pick k (or about k) entries uniformly at random
[I, J] = find(W);
% [I, J, k] = randmask(m, n, k);

% Compute the values of AB at these entries (this is a C-Mex function)
% M(W==0)=0;
X = M(W==1);
% X = spmaskmult(A, B, I, J);

% Define the confidence we have in each measurement X(i) IMPORTANT: If you
% wish to use non-uniform weights (so, C is not a constant vector), you
% need the latest compiled version of spbuildmatrix.c. You can obtain it by
% executing:
%   mex -lmwlapack -lmwblas -largeArrayDims spbuildmatrix.c
% You do not need to do that if you are using .mexw64 or .mexmaci64 files
% (they are already up to date). If you compile the file for a new
% platform, it would be great if you sent us the output of that
% compilation, for inclusion in the package: nicolasboumal@gmail.com.
% Thanks!
C = ones(size(X));

% Add noise if desired noisestd = 0; X = X + noisestd*randn(size(X));


%% Feeding the problem instance to RTRMC

% Pick a value for lambda, the regularization parameter
% lambda = 0;


% perm = randperm(k); I = I(perm); J = J(perm); X = X(perm); C = C(perm);


% Build a problem structure
problem = xt_nb_buildproblem(I, J, X, C, m, n, r, opts.nu);


% Compute an initial guess
% initstart = tic;
% U0 = initialguess(problem);
% U0 = je_generate_sample(sample_num, m, n, r);
% inittime = toc(initstart);

% [Optional] If we want to track the evolution of the RMSE as RTRMC
% iterates, we can do so by specifying the exact solution in factored form
% in the problem structure and asking RTRMC to compute the RMSE in the
% options structure. See the subfunction computeRMSE in rtrmc.m to see how
% to compute the RMSE on a test set if the whole matrix is not available in
% a factorized A*B form (which is typical in actual applications).

% Setup the options for the RTRMC algorithm These are the algorithms shown
% in the papers:
%  RTRMC 2  : method = 'rtr', order = 2, precon = false RTRMC 2p : method =
%  'rtr', order = 2, precon = true RTRMC 1  : method = 'rtr', order = 1,
%  precon = false RCGMC    : method = 'cg',  order = 1, precon = false
%  RCGMCp   : method = 'cg',  order = 1, precon = true
opts.method = 'rtr';     % 'cg' or 'rtr', to choose the optimization algorithm
opts.order = 2;          % for rtr only: 2 if Hessian can be used, 1 otherwise
opts.precon = true;      % with or without preconditioner
opts.maxiter = opts.max_iter;      % stopping criterion on the number of iterations
opts.maxinner = opts.max_trials;      % for rtr only : maximum number of inner iterations
opts.tolcost = opts.tol; % stopping criterion on the cost function value
% opts.tolgradnorm = opts.tol; % stopping criterion on the norm of the gradient
opts.verbosity = opts.display * 2; % how much information to display during iterations
opts.computeRMSE = true; % set to true if RMSE is to be computed at each step


% Call the algorithm here. The outputs are U and W such that U is
% orthonormal and the product U*W is the matrix estimation. The third
% output, stats, is a structure array containing lots of information about
% the iterations made by the optimization algorithm. In particular, timing
% information should be retrieved from stats, as in the code below.
[U, V, stats] = xt_nb_rtrmc(problem, opts, U);
V = V';

% time = inittime + [stats.time];
% rmse = [stats.RMSE];
iter = [stats.iter];
iter = iter(end);

cost = norm(W .* (M - U * V'), 'fro')^2;
% switch opts.method
%     case 'rtr'
%         semilogy(time, rmse, '.-');
%     case 'cg'
%         semilogy(time, rmse, '.-'); S = [stats.linesearch]; ls_evals =
%         [S.costevals]; for i = 1 : length(ls_evals)
%             text(time(i), rmse(i)*1.20, num2str(ls_evals(i)));
%         end title(sprintf('Numbers indicate line search cost
%         evaluations.'));
% end xlabel('Time [s]'); ylabel('RMSE');

end
