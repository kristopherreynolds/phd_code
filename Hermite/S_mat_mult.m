function [G,S,E] = S_mat_mult(N,t)
%function [G,S,E] = S_mat_mult(N,t)
%

S0 = [];
load('S45_500N.mat');

S = cell(N+1,1);
G = cell(N+1,1);
E = cell(N+1,1);

hwb = waitbar(0,'Doing Matrix Multiplication for S');
for m = 0 : N
    Ecurrent = diag(1i.^(0:m));
    Gcurrent = diag(exp(1i.*(m:-2:-m)*t));
    S0current = S0{m+1};
    Scurrent = real(Ecurrent*S0current'*Gcurrent*S0current*conj(Ecurrent));
    S{m+1} = Scurrent;
    G{m+1} = Gcurrent;
    E{m+1} = Ecurrent;
    waitbar(m/N)
end
close(hwb)