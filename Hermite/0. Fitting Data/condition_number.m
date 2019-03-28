function out = condition_number(R,N,npts)

inds = 1 : npts;
x = -R+2*R*(inds-1)/(npts-1);
H = herm_funs(N,x);
[~,S,~] = svd(H);
s = diag(S);
s = s(1:N+1);
c = max(s)/min(s);

out = c;