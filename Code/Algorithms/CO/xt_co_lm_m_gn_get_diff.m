function [b, H, N_perp]=xt_co_lm_m_gn_get_diff(M, mask, N, r)
%N=randn(9,3);M0=randn(20,3)*N';M=M0+randn(20,9)*0.1;r=3;mask=rand(20,9)>0.3;N=N+0.01*randn(9,3);N=orth(N);
[m,n]=size(M);
b=zeros((n-r)*r,1);
H=zeros((n-r)*r);
% b=zeros(n*r,1);
% H=zeros(n*r);
N_perp=null(N');
for i=1:m
  if(mod(i,100)==0)
    j=0;
  end
  mm=mask(i,:)';
  index=find(mm==1);
  n0=length(index);
  I=eye(n);
  L=N_perp(index,:);
  %     L=I(index,:);
  D=kron(eye(r),L);
  
  mi=M(i,index)';
  Ni=N(index,:);
  E=inv(Ni'*Ni);
  A=Ni*E;
  P=A*Ni';
  C=xt_co_communication(n0,r);
  
  temp1=A'*mi;
  temp2=(eye(n0)-P)*mi;
  
  fi=temp2;
  Bi=(kron(temp1',eye(n0)-P)+kron(temp2',A)*C)*D;
  
  H=H+Bi'*Bi;
  b=b-Bi'*fi;
  
  
  %     b=b-D'*(kron(temp1,temp2)+C'*kron(temp2,temp1));
  %
  % %     temp=C'*kron(temp2,eye(r))*(kron(temp2',E)*C-2*kron(temp1',A'));
  % %     temp=temp-kron(temp1,eye(n0))*kron(temp1',eye(n0)-P);
  % %     temp=temp+temp';
  %
  %     temp=kron(eye(r),temp2)*(kron(temp2',E)*C-2*kron(temp1',A'));
  %     temp=temp-C'*kron(eye(n0),temp1)*kron(temp1',eye(n0)-P);
  %     temp=temp+temp';
  %     H=H-D'*temp*D;
end
H=H;
b=-b;

end