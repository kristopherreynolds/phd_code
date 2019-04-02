function out = thi_num(m,n,k,a,b,c,x)
%out = thi_num(m,n,k,a,b,c,x)

hm = herm_fun(m,a*x);
hn = herm_fun(n,b*x);
hk = herm_fun(k,c*x);

dx = mean(diff(x));

out = sum(hm.*hn.*hk.*dx);