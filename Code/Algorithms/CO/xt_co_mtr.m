function M = xt_co_mtr(v, m, n)
M=[];
[k,l]=size(v);
if(k==1)
  v=v';
end
if(length(v)==m*n)
  for i=1:n
    M=[M,v(1+m*(i-1):m*i)];
  end
end