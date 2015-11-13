function [e, no, N] = xt_co_lm_m(M, mask, r, N, max_iter, tol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LM_M in the paper
%  M \in  R^{m,n} with m>=n. mask with 1 for observed 0 for missing.
%  Entries in M can't be 'NaN', should be any real value, 0, or random
%  value. r is the rank. N is the initialization. N \in R^{n,r}, i.e., the
%  initial subspace of row vectors
% Output the objectives at each iteration. no is the iteration number.
% Ns are the estimated subspaces at each iteration.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 5, max_iter = 300;
end
if nargin < 6, tol = 1e-10;
end

%M=randn(20,9);r=3;mask=rand(20,9)>0.3;N=rand(9,3);%N=orth(N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% When used, pls note this part.
% N=randn(9,3);M0=randn(20,3)*N';M=M0+randn(20,9)*0.2;r=3;mask=rand(20,9)>0.3;N=N+0.01*randn(9,3);N=orth(N);
% [m,n]=size(M);
% if(size(N)==[m,n])
%     [u,s,v]=svds(N,r);
% else
%     if(size(N)==[n,r])
%         N=orth(N);
%     else
%         return;
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for simulation
% N=randn(9,3);M0=randn(20,3)*N';M=M0+randn(20,9)*0.2;r=3;mask=rand(20,9)>0.3;N=N+0.01*randn(9,3);N=orth(N);
% [m,n]=size(M);
% if(size(N)==[m,n])
%     [u,s,v]=svds(N,r);
% else
%     if(size(N)==[n,r])
%         N=orth(N);
%     else
%         return;
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[m,n]=size(M);

converged=0;
times=0;
e(1)=xt_co_obj_missing(M,mask,r,N);
e_old=e(1);
lambda=0.001;
% lambda=1e-6;
k=10;
N_old=N;
while converged==0
  [b,H,N_com]=xt_co_lm_m_get_diff(M,mask,N_old,r);
  
  while(1)
    delta=inv(H+lambda*eye(size(H,1)))*b;
    %delta=inv(H+lambda*diag(diag(H)))*b;
    
    % variation=b'*delta;
    K=xt_co_mtr(delta,n-r,r);
    N_new=N_old+N_com*K;
    e_new=xt_co_obj_missing(M,mask,r,N_new);
    
    if(e_new<e_old || abs(e_new-e_old)<tol)
      lambda=lambda/k;
      break;
    else
      lambda=lambda*k;
    end
    
  end
  
  N_new=orth(N_new);
  d=svd(N_old'*N_new);
  %if(e_old-e_new>=0&e_old-e_new<1e-10)
  if(1-min(d)<1e-10)
    converged=1;
  end
  
%   if(mod(times,10)==0)
%     tmp=times/10+1;
%     Ns(:,:,tmp)=N_new;
%   end
  
  %     lambda
  %     e_old-e_new
  e(times+2)=e_new;
  N_old=N_new;
  e_old=e_new;
  times=times+1;
  if(times==max_iter)
    converged=1;
  end
end
%e(1)=times;
no=times;
i=0;
% Ns(:,:,tmp+1)=N_new;
N=N_new;
% times
% e_new
% if(sqrt(min(e)/sum(sum(mask)))<0.0226)
%     rms=sqrt(min(e)/sum(sum(mask)));
%     save('N','N','rms');
% end

end