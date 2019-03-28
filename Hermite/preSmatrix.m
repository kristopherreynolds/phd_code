function out = preSmatrix(q,n,m)
%function out = preSmatrix(q,n,m)

kd = @(m,n) m==n;
out = sqrt(n)*sqrt(m-n+1)*kd(q,n-1) - sqrt(n+1)*sqrt(m-n)*kd(q,n+1);