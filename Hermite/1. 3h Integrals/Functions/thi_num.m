function out = thi_num(m,n,k,a,b,c,x)
%out = thi_num(m,n,k,a,b,c,x)


if  round((m+n+k)/2) == (m+n+k)/2
hm = herm_fun(m,a*x);
hn = herm_fun(n,b*x);
hk = herm_fun(k,c*x);
dx = mean(diff(x));
out = sum(hm.*hn.*hk.*dx);
else    
out = 0;    
end