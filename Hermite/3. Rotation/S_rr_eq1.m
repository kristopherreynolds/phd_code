function out = S_rr_eq1(S,n,q,m,t)
%function out = S_rr_eq1(S,n,q,m,t)
%This recurrence relation gives S_{n,q}^{m+1} and S_{n+1,q}^{m+1}
%from S_{n,q}^m and S_{n,q-1}^m
%The input S is a 1x(m+1) cell

if q-1 < 0
    out = 0;
else
    %inputs
    nqm = S{m+1}(n+1,q+1);
    nq_m1m = S{m+1}(n+1,q-1+1);
    
    %outputs
    A1 = [1/sqrt(n+1) 0; 0 1/sqrt(m-n+1)];
    B1 = [cos(t) -sin(t);sin(t) cos(t)];
    
    ksi_p1 = A1*B1*[sqrt(q)*nq_m1m;sqrt(m-q+1)*nqm];
    
    out.svalues = ksi_p1;
    out.indices = [n+1 q m+1; n q m+1]+1;
end
