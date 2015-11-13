function f = xt_co_obj_missing(M, mask, r, N)

%N=randn(9,3);M0=randn(20,3)*N';M=M0+randn(20,9)*0.1;r=3;mask=rand(20,9)>0.3;N=N+0.01*randn(9,3);N=orth(N);
[m,n]=size(M);

f=0;
for i=1:m
  mm=mask(i,:)';
  index=find(mm==1);
  
  mi=M(i,index)';
  Ni=N(index,:);
  
  %     P=Ni*inv(Ni'*Ni)*Ni';
  %     temp=mi-P*mi;
  %     f=f+mi'*temp;
  b=Ni\mi;
  temp=mi-Ni*b;
  f=f+temp'*temp;
end

end