function C = communication(p, q)
% vec(A')=C*vec(A)  with A in R^{p,q}
temp=0:p:p*(q-1);
temp=temp';
index=[];
for i=1:p
  index=[index;temp+i];
end
C=eye(p*q);
C=C(index,:);

end