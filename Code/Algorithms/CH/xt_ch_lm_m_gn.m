function [e, iter, Ns] = xt_ch_lm_m_gn(M, mask, r, N, max_iter, tol, disp_level)

% LM_M_GN (Damped RW1 + ambiguity constraint A)

if nargin < 5, max_iter = 300;
end
if nargin < 6, tol = 1e-10;
end
if nargin < 7, disp_level = 1;
end

vec = @(X) X(:);
[m,n]=size(M);
n_r = n - r;
nr_r2 = n_r * r;

% Form cells of repeatedly-used variables.
store.C = cell(m, 1);
store.index = cell(m, 1);
store.n0 = cell(m, 1);
store.rng = cumsum([0; sum(mask, 2)]);
for i = 1 : m
  store.index{i} = find(mask(i,:)' == 1);
  store.n0{i} = length(store.index{i});
  store.C{i} = sparse(1 : (store.n0{i} * r), vec(reshape(1 : (store.n0{i} * r), store.n0{i}, r)'), 1);
  % store.C{i} = communication(store.n0{i}, r);
end

% converged=0;
% times=0;
[N, ~] = qr(N, 0);
[e, NiQ, NiR] = obj_missing(M, N, store);
lambda=1e-4;
% lambda=1e-6;
% k=10;
k_inc = 10;
k_dec = 100;
max_trials = 150;

for iter = 1 : max_iter
  [b, H, N_perp] = get_b_H(M, N, NiQ, NiR, r, store);
  
  for trial = 1 : max_trials
    delta = (H + lambda * speye(nr_r2)) \ b;
    %delta=inv(H+lambda*diag(diag(H)))*b;
    
    % variation=b'*delta;
    % K=mtr(delta,n-r,r);
    [N_new, ~] = qr(N + N_perp * reshape(delta, n_r, r), 0);
    [e_new, NiQ, NiR] = obj_missing(M, N_new, store);
    
    if e_new < e
      lambda = lambda / k_dec;
      break;
    else
      lambda = lambda * k_inc;
    end
  end

  if disp_level
    fprintf('[%04d] %.2e %.4e\n', iter, lambda, e_new);
  end
  
  if (trial == max_trials) || (abs(e - e_new) < tol * e), break
  end
  
  N = N_new;
  e = e_new;
  
end

Ns = N_new;
e = e_new;

end

function [b, H, N_perp] = get_b_H(M, N, NiQ, NiR, r, store)

[m,n]=size(M);
b=zeros(n*r, 1);
H=zeros(n*r);
N_perp=null(N');

for i = 1 : m
  % if(mod(i,100)==0)
  %   j=0;
  % end
  % mm=mask(i,:)';
  % index=find(mm==1);
  % n0=length(index);
  C = store.C{i};
  index = store.index{i};
  n0 = store.n0{i};
  
  % I=eye(n);
  % L=N_perp(index,:);
  %     L=I(index,:);
  % D=kron(eye(r),L);
  rng = bsxfun(@plus, index - 1, 1:n:n*r);
  rng = rng(:);
  
  mi=M(i,index)';
  % Ni=N(index,:);
  % E=inv(Ni'*Ni);
  % A=Ni*E;
  % P=A*Ni';
  % C=communication(n0,r);
  
  E = eye(r) / NiR{i};
  % AT = NiR{i} \ NiQ{i}';
  
  % temp1=A'*mi;
  % temp2 = (eye(n0)-P)*mi;
  temp0 = NiQ{i}' * mi;
  temp1 = NiR{i} \ temp0;
  temp2 = mi - NiQ{i} * temp0;
  b(rng) = b(rng) + kron(temp1, temp2);
  
  % temp = kron(temp1 * temp1', speye(n0) - NiQ * NiQ'); % RW2
  % temp = temp + C' * kron(temp2 * temp2', E * E') * C; % RW1
  % temp = temp + C' * kron(temp2,eye(r)) * (2 * kron(temp1',AT) - kron(temp2',E) * C); % RW0
  % temp = C' * kron(temp2 * temp1', AT);
  H(rng, rng) = H(rng, rng) ...
    + kron(temp1 * temp1', speye(n0) - NiQ{i} * NiQ{i}') ...
    + C' * kron(temp2 * temp2', E * E') * C;
end

PT = kron(speye(r), N_perp');
b = PT * b;
H = PT * H * PT';
H = 0.5 * (H + H');

end

function [f, NiQ, NiR] = obj_missing(M, N, store)

m = size(M, 1);
NiQ = cell(m, 1);
NiR = cell(m, 1);

f=0;
for i=1:m
  % mm=mask(i,:)';
  % index=find(mm==1);
  index = store.index{i};
  
  mi=M(i,index)';
  % Ni=N(index,:);
  [NiQ{i}, NiR{i}] = qr(N(index, :), 0);
  
  %     P=Ni*inv(Ni'*Ni)*Ni';
  %     temp=mi-P*mi;
  %     f=f+mi'*temp;
  % b=Ni\mi;
  % temp=mi-Ni*b;
  temp = mi - NiQ{i} * (NiQ{i}' * mi);
  f=f+temp'*temp;
end

end
