function [Sout] = S_mat_mult(N,t)
%function [S] = S_mat_mult(N,t)
%

S0 = [];
load('steering_coeffs_45deg.mat');
S45 = S0;
Sout = cell(N+1,1);

hwb = waitbar(0,'Doing Matrix Multiplication for S');
for m = 0 : N
    E = diag(1i.^(0:m));
    G = diag(exp(1i*(m:-2:-m)*t));
    S0 = S45{m+1};
    Sout{m+1} = real(E*S0'*G*S0*E');
    waitbar(m/N)
end
close(hwb)