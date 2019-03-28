function out = S_rr_eq2(S,n,q,m,t)
%function out = S_rr_eq2(S,n,q,m,t)
%This recurrence relation gives S_{n,q}^{m+1} and S_{n,q+1}^{m+1}
%from S_{n,q}^m and S_{n-1,q}^m
%The input S is a 1x(m+1) cell

if n-1 < 0
    out = 0;
else
    %inputs
    nqm = S{m+1}(n+1,q+1);
    n_m1qm = S{m+1}(n-1+1,q+1);
    
    %outputs
    A2 = [1/sqrt(q+1) 0; 0 1/sqrt(m-q+1)];
    B2 = [cos(t) sin(t);-sin(t) cos(t)];
    
    ksi_p1 = A2*B2*[sqrt(n)*n_m1qm;sqrt(m-n+1)*nqm];
    
    out.svalues = ksi_p1;
    out.indices = [n q+1 m+1; n q m+1]+1;
end
