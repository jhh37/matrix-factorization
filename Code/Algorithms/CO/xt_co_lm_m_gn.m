function [e, no, N] = xt_co_lm_m_gn(M, mask, r, N, max_iter, tol)

if nargin < 5, max_iter = 300;
end
if nargin < 6, tol = 1e-10;
end

%M=randn(20,9);r=3;mask=rand(20,9)>0.3;N=rand(9,3);%N=orth(N);
%N=randn(9,3);M0=randn(20,3)*N';M=M0+randn(20,9)*0.2;r=3;mask=rand(20,9)>0.3;N=N+0.01*randn(9,3);N=orth(N);
%[m,n]=size(M);
% if(size(N)==[m,n])
%     [u,s,v]=svds(N,r);
% else
%     if(size(N)==[n,r])
%         N=orth(N);
%     else
%         return;
%     end
% end

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
  %     [b,H,N_com]=get_b_H(M,mask,N_old,r);
  [b,H,N_com]=xt_co_lm_m_gn_get_diff(M,mask,N_old,r);
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
% function u = CLM (X1, Y1, X2, Y2, F0, u0)
%
% 	 iter = 0;
% 	 iter_max = 100;
% 	 c = 0.0001;
%          eps = 1.0e-15;
%
% 	 Xi = calcXi (X1, Y1, X2, Y2, F0);
%
% 	 F = reshape (u0, 3, 3);
% 	 F = F';
% 	 [U, S, V] = svd (F);
% 	 theta = asin (S(2,2) / sqrt (S(1,1) * S(1,1) + S(2,2) * S(2,2)));
% 	 S(1, 1) = cos (theta);
% 	 S(2, 2) = sin (theta);
% 	 S(3, 3) = 0.0;
% 	 F  = U * S * V';
% 	 F_ = eye (3, 3);
% 	 J = calcJ (F, Xi, X1, Y1, X2, Y2, F0);
%
% 	 while (iter < iter_max)
%
% 	       iter = iter + 1;
%
% 	       Fu = calcFu (F);
% 	       Fv = calcFv (F);
% 	       ut = calcUt (U, S, V);
% 	       u  = reshape (F', 9, 1);
%
% 	       X = calcX (u, Xi, X1, Y1, X2, Y2, F0);
%
% 	       [H, g] = setHandG (X, Fu, Fv, u, ut);
% 	       DH = diag (diag(H), 0);
%
% 	       while (1)
% 		     param = (H + c * DH) \ (-g);
% 		     omega = param(1:3, 1);
% 		     U_ = createRotationMatrix (omega) * U;
% 		     omega = param(4:6, 1);
% 		     V_ = createRotationMatrix (omega) * V;
% 		     theta_ = theta + param(7);
% 		     S_ = zeros(3, 3);
% 		     S_(1, 1) = cos (theta_);
%      		     S_(2, 2) = sin (theta_);
% 		     F_ = U_ * S_ * V_';
% 		     J_ = calcJ (F_, Xi, X1, Y1, X2, Y2, F0);
%
% 		     if (J_ < J | abs (J_ - J) < eps)
%                  J = J_;
%                  c = c * 0.1;
%                  break;
%              else
%                  c = c * 10.0;
%              end
%            end
%
% 	       if (calcDistance (F, F_) < 0.5e-12)
%                u = reshape (F_', 9, 1);
%                break;
%            else
%                F = F_;
%                U = U_;
%                S = S_;
%                V = V_;
%                theta = theta_;
%            end
% 	 end
% endfunction

end