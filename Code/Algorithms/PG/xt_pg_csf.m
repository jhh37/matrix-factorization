function [M, S, X, rmse, iter] = xt_pg_csf(W, r, X0, numIter, RMSE_TOL, verbose, B)
%function [M,S,X,rmse,iter] = pgCSF ( W, r, X0, numIter, RMSE_TOL, verbose, B )
%
% By Paulo Gotardo
%
% Computes factorization W = MS, with M = BX
%
% Inputs:
%
% W is the observation matrix (missing data encoded as NaN entries)
% r is the factorization rank
% X0 is the initial value of X (or the initial M0 when basis B is the identity)
% numIter is the maximum number of iterations
% RMSE_TOL is the parameter of the stopping (convergence) rule
% verbose is boolean
% B is the basis of M = BX (default basis is the identity matrix)
%
xt_pg_check_buildkron()

if (nargin < 7), B = eye(size(W,1)); end  % assume canonical basis I_m

% Auxiliary variables
VALID = isfinite(W);                 % mask of observed data
% (missing data marked as NaN in W)
rmse  = zeros(numIter+1,1);          % error at each iteration

n = size(W,2);                       % number of points
d = size(B,2);                       % number of DCT basis vectors
nx  = d*r;                           % number of unknowns
Idr = speye(nx);                     % damping matrix
delta = 1e-4;                        % initial damping parameter
g   = zeros(nx,1);                   % gradient vector
H   = zeros(nx,nx);                  % Hessian matrix

% Initialize X
if isempty(X0)
  X = [ eye(r,r) ; zeros(d-r,r) ]; % deterministic initialization
else
  [X, ~] = qr(X0, 0);
  if (size(X,1) < d), X(d,:) = 0; end  % enlarge
end

warning('off', 'MATLAB:nearlySingularMatrix');     % OK! Damping fixes Hessian

% Compute initial factors
M = B*X;                             % current estimate of M
S = zeros(r,n);                      % current estimate of S
R = zeros(size(W));                  % residual values
for j = 1:n
  mask = VALID(:,j);
  wj = W(mask,j);
  Mj = M(mask,:);
  %S(:,j) = pinv(Mj) * wj;
  S(:,j) = Mj \ wj;
  R(mask,j) = wj - Mj * S(:,j);
end

% Compute initial fit error (initial cost f(X))
rmse(1) = sqrt( nanmean( R(VALID(:)).^2 ));
if (verbose)
  fprintf('\n\ni = 0 \t RMSE = %-15.10f \n', rmse(1) )
end

% Main loop

for iter = 1:numIter
  
  % (1) calculate Gradient and Jacobian (J'J) approx to Hessian (Gauss-Newton)
  g(:) = 0; H(:) = 0;
  for j = 1:n
    mask = VALID(:,j);
    Mj = M(mask,:);
    Bj = B(mask,:);
    sj = S(:,j);
    rj = R(mask,j);
    
    %PjBj = Bj - Mj * (pinv(Mj)*Bj);            % (I-Pj)*Bj
    PjBj = Bj - Mj * (Mj\Bj);                  % (I-Pj)*Bj
    Jj = xt_pg_kronmex(sj', PjBj);
    g = g - (rj' * Jj)';
    %H = H + Jj'*Jj;                            % slow
    H = H + xt_pg_kronmex(sj*sj', Bj'*PjBj);         % faster
    
    %g = g - reshape((Bj'*rj)*sj', [], 1);
  end
  H = 0.5 * (H + H');                            % force H symmetric
  
  % (2) Repeat solving for vec_dX until f(X-dX) < f(X) or converged
  while true
    %vec_dX = pinv(H + delta*Idr) * g;
    vec_dX = (H + delta*Idr) \ g;
    newX = X - reshape(vec_dX, d, r);
    % [newX,foo] = qr(newX,0);
    [newX, ~] = qr(newX,0);
    
    % Compute new factors
    M = B * newX;
    for j = 1:n
      mask = VALID(:,j);
      wj = W(mask,j);
      Mj = M(mask,:);
      %S(:,j) = pinv(Mj) * wj;
      S(:,j) = Mj \ wj;
      R(mask,j) = wj - Mj * S(:,j);
    end
    
    % Evaluate cost f(newX)
    R2 = R(VALID(:)).^2;
    max_err = sqrt(nanmax( R2(:) ));
    rmse(iter+1) = sqrt(nanmean( R2(:) ));
    
    % Damping termination tests
    if (rmse(iter+1) < rmse(iter)), OK = true; break, end
    
    % Continue damping of H?
    delta = delta * 10;
    if (delta > 1.0e30), OK = false; break, end
  end                                                            % end while
  % Error test (bailed out with no descent)
  if (~OK), disp('Error: cannot find descent direction!'), break, end
  
  % (3) Book-keeping; display new errors
  delta = max( delta / 100, 1.0e-20);
  X = newX;
  
  % (4) Display new error and test convergence
  if (verbose)
    fprintf('i = %-4d  RMSE = %-9.6f (max %-9.6f)  l = 1.0e%03d\n', ...
      iter, rmse(iter+1), max_err, fix(log10(delta)) )
  end
  if (rmse(iter) - rmse(iter+1) < RMSE_TOL), break, end
end

% Truncate vector of RMSE values
iter = iter + 1;
rmse = rmse(iter);

warning('on', 'MATLAB:nearlySingularMatrix');     % OK! Damping fixes Hessian

% ----------------------------------------------------------------------------

end