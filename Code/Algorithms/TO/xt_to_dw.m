function [U, V, err, iter, err_log, lambda_log] = xt_to_dw(Y, H, r, Vini, tol, max_iter, disp_level)
% The damped Wiberg algorithm for low-rank matrix factorization Y -> U V'
%
% [U, V] = damped_wiberg(Y, H, r) factorizes the data matrix Y (say m x n)
% into the product of U (m x r) and V (n x r). Note that r specifies the
% column size of U and V. The indicator H that has the same size as Y
% specifies non-missing/missing components, i.e., H(i, j) = 0 if Y(i, j) is
% missing and H(i, j) = 1 if Y(i, j) is non-missing.
%
% [U, V] = damped_wiberg(Y, H, r, Vini) uses Vini as initial values of V
% instead of random initial values.
%
% [U, V] = damped_wiberg(Y, H, r, Vini, tol, max_iter) further uses tol
% instead of the default torelance (1e-6) and max_iter instead of the
% default maximum iteration counts (500).
%
% The implementation is based on:
% [1] Takayuki Okatani, Takahiro Yoshida, Koichiro Deguchi: Efficient
%     algorithm for low-rank matrix factorization with missing components
%     and performance comparison of latest algorithms. Proc, ICCV 2011:
%     842-849.

if nargin < 4, Vini = randn(size(Y,2), r);
end
if nargin < 6, max_iter = 100;
end
if nargin < 5, tol = 1e-6;
end
if nargin < 7, disp_level = 1;
end

%==========================================================================
% Check dimensions
%==========================================================================
%m = size(Y,1);
%n = size(Y,2);
[m, n] = size(Y);

% Make the vector u & v from the matrix U & V
Vt = Vini.';
v = Vt(:);

%==========================================================================
% Count the number of nonmissing components
%==========================================================================
p = 0;
Findx = zeros(m+1, 1);
count = 1;
for i=1:m
  Findx(i) = count;
  for j=1:n
    if H(i,j) > 0,
      count = count + 1;
      p = p + 1;
    end
  end
end
Findx(m+1) = p+1;

%disp(sprintf('Number of nonmissing components = %d\n', p));

% F = zeros(p, m*r);
% G = zeros(p, n*r);
y = zeros(p, 1);

% F -> FQ*FR (QR-decomposition)
Fpack = zeros(p, r);    % dense matrix
FQpack = zeros(p, r);   % dense matrix
FRpack = zeros(m*r, r);  % submatrices are upper-triangular

% Make y from Y
q = 1;
for i=1:m
  for j=1:n
    if H(i,j) > 0.5,
      y(q) = Y(i,j);
      q = q+1;
    end
  end
end

lambda = 0.1;

%log
err_log = zeros(max_iter,1);
lambda_log = zeros(max_iter,1);

%==========================================================================
% Start of the loop
%==========================================================================
err_last = 0.;
for iter = 1:max_iter
  
  % orthogonalize v
  
  % if 0
  %   V1 = zeros(n, r);
  %   for j=1:n
  %     V1(j,:) = v((j-1)*r+1:j*r)';
  %   end
  %   [V1, V1r] = qr(V1,0);
  %   %    v = V1'(:);
  %   for j=1:n
  %     v((j-1)*r+1:j*r) = V1(j,:)';
  %   end
  % end
  
  q = 1;
  for i=1:m
    for j=1:n
      if H(i,j) > .5,
        Fpack(q, :) = v((j-1)*r+1:(j-1)*r+r);
        q = q + 1;
      end
    end
  end
  
  % QR decomposition of F
  for i=1:m
    [FQpack(Findx(i):Findx(i+1)-1,1:r), FRpack((i-1)*r+1:(i-1)*r+r,1:r)] = qr(Fpack(Findx(i):Findx(i+1)-1,1:r), 0);
  end
  
  % 2.2 Compute u minimizing |Fu - y|^2
  [u, err, errvec] = update_u(Fpack, Findx, FQpack, FRpack, y, m, r);
  
  %disp(sprintf('norm u= %f norm v=%f', norm(u), norm(v)));
  
  % 3. Check convergence
  if disp_level ~= 0
    err_log(iter) = err;
    lambda_log(iter) =  lambda;
    fprintf('[%04d] %.2e %.4e\n', iter, lambda, err);
    drawnow
  end
  if abs(err - err_last) < tol*err, break
  end
  err_last = err;
  
  % Make G'QFG = G'G - (Q1'G)'(Q1'G)
  Gyvec = zeros(n*r,1);
  % Gyvec1 = zeros(m*r,1);
  %GQG = zeros(n*r,n*r);
  QG  = zeros(m*r,n*r);%sparse(m*r,n*r);
  
  
  for i=1:m
    q =1;
    FQpack2 = FQpack(Findx(i):Findx(i+1)-1,:)';
    FQv = FQpack2(:);
    ui = u((i-1)*r+1:i*r);
    tmp2 = FQv * ui';
    for j=1:n
      if H(i,j) ~= 0,
        %                 tmp2 = FQpack(Findx(i)-1+q,:)'*u((i-1)*r+1:i*r)';
        %                 for l=1:r
        %                     for k=1:r
        %                         QG((i-1)*r+k,(j-1)*r+l) = tmp2(k,l);
        %                     end
        %                 end
        %QG((i-1)*r+1:i*r,(j-1)*r+1:j*r) = FQpack(Findx(i)-1+q,:)'*u((i-1)*r+1:i*r)';
        QG((i-1)*r+1:i*r,(j-1)*r+1:j*r) = tmp2((q-1)*r+1:q*r,:);
        %Gyvec((j-1)*r+1:j*r) = Gyvec((j-1)*r+1:j*r) + y(Findx(i)-1+q) * ui;
        % Gyvec = G'*e;
        Gyvec((j-1)*r+1:j*r) = Gyvec((j-1)*r+1:j*r) + errvec(Findx(i)-1+q) * ui;
        q = q + 1;
      end
    end
    %Gyvec1((i-1)*r+1:i*r) = FQpack(Findx(i):Findx(i+1)-1,1:r)'*y(Findx(i):Findx(i+1)-1);
    %Gyvec1((i-1)*r+1:i*r) = FQpack2*y(Findx(i):Findx(i+1)-1);
  end
  
  %    GQG = -QG'*QG;
  %    Gyvec = Gyvec-QG'*Gyvec1;
  QGs = sparse(QG);
  GQG = -QGs'*QGs;
  %Gyvec = Gyvec-QGs'*Gyvec1;
  for j=1:n
    tmp1 = zeros(r,r);
    %q = 1;
    for i=1:m
      if H(i,j) > 0.5,
        tmp1 = tmp1 + u((i-1)*r+1:i*r)*u((i-1)*r+1:i*r)';
      end
    end
    GQG((j-1)*r+1:j*r,(j-1)*r+1:j*r) = GQG((j-1)*r+1:j*r,(j-1)*r+1:j*r) + tmp1;
  end
  
  % Make M from v
  M = zeros(n*r, r*r);
  for i=1:r
    for j=1:r
      for k=1:n
        M((k-1)*r+j, (i-1)*r+j) = v((k-1)*r+i);
      end
    end
  end
  [Mq, ~] = qr(M,0);
  
  % 4.2 Determine dv minimizing |QF G dv - QF y|^2
  
  gauss_newton = 0;
  while 1
    % Amat = zeros(size(G,2),size(G,2));
    if gauss_newton == 1,
      Amat = GQG + Mq*Mq';
    else
      % Marquardt version
      %Amat = GQG + Mq*Mq' + lambda*diag(diag(GQG + Mq*Mq'));
      % original (Levenberg)
      Amat = GQG + Mq*Mq' + lambda*diag(ones(size(Mq,1),1));
      
      %Amat = G'*QFG + Mq*Mq' + lambda*diag(ones(size(Mq,1),1));
      %Amat = GQG + M*M' + lambda*diag(ones(size(Mq,1),1));
      %Amat = GQG + lambda*diag(ones(size(Mq,1),1));
      %Amat = GQG + Mq*Mq' + lambda*(diag(ones(size(Mq,1),1)) - Mq*Mq');
      %bvec = zeros(size(Amat,1),1);
      %bvec = G'*QFy;
    end
    bvec = Gyvec;
    
    %{
        R = chol(Amat);
        dv3 = R'\bvec;
        dv2 = R\dv3;
    %}
    dv2 = Amat \ bvec;
    
    % v1 = zeros(size(v));
    v1 = v + dv2;
    
    % 2.1 Make F from v
    q = 1;
    for i=1:m
      for j=1:n
        if H(i,j) > 0,
          Fpack(q, :) = v1((j-1)*r+1:(j-1)*r+r);
          q = q + 1;
        end
      end
    end
    
    % QR decomposition of F
    for i=1:m
      [FQpack(Findx(i):Findx(i+1)-1,1:r), FRpack((i-1)*r+1:(i-1)*r+r,1:r)] = qr(Fpack(Findx(i):Findx(i+1)-1,1:r), 0);
    end
    
    % 2.2 Compute u minimizing |Fu - y|^2
    [u, err] = update_u(Fpack, Findx, FQpack, FRpack, y, m, r);
    
    if gauss_newton == 1
      break;
    end
    
    %disp(sprintf('[***] %e %e', lambda, err));
    
    if err > err_last,
      lambda = lambda * 10;
    else
      lambda = lambda * 0.1;
      break;
    end
  end
  
  %disp(sprintf('lambda = %e', lambda));
  
  % 4.3 Update v
  v = v + dv2;
end

% Recover U and V from u and v
U = zeros(m, r);
V = zeros(n, r);
for i=1:m
  U(i,:) = u((i-1)*r+1:(i-1)*r+r);
end
for j=1:n
  V(j,:) = v((j-1)*r+1:(j-1)*r+r);
end


function [u, err, errvec] = update_u(Fpack, Findx, Q1pack, R1pack, y, m, r)
% Solve |Fu-y|^2 -> min for u
% R1 u = Q1' y
% F = F(v) = (expansion of Q1pack)*(exp. of R1pack)

u = zeros(m*r,1);

for i=1:m
  %disp(sprintf('%e', R1pack((i-1)*r+1,1)/R1pack((i-1)*r+r,r)));
  
  x = R1pack((i-1)*r+1:i*r,1:r) \ Q1pack(Findx(i):Findx(i+1)-1,1:r)'*y(Findx(i):Findx(i+1)-1);
  u((i-1)*r+1:i*r) = x;
end

errvec = zeros(size(y));
for i=1:m
  errvec(Findx(i):Findx(i+1)-1) = y(Findx(i):Findx(i+1)-1) - Fpack(Findx(i):Findx(i+1)-1,1:r)*u((i-1)*r+1:i*r);
end
err = errvec.'*errvec;
