function [G,S0,S,E] = S_mat_mult(N,t)
%function [G,S0,S,E] = S_mat_mult(N,t)
%

if N<=100
S0 = load('S45_100N.mat');
S0 = S0.S_45_100N;
elseif N<=200
S0 = load('S45_200N.mat');
S0 = S0.S_45_200N;
elseif N<=300
S0 = load('S45_300N.mat');
S0 = S0.S_45_300N;
elseif N<=400
S0 = load('S45_400N.mat');
S0 = S0.S_45_400N;
elseif N<=500
S0 = load('S45_500N.mat');
S0 = S0.S_45_500N;
elseif N<=600
S0 = load('S45_600N.mat');
S0 = S0.S_45_600N;
else
error('Exceeds Bandlimit of pre-stored S Matrix')
end


S = cell(N+1,numel(t));
G = cell(N+1,numel(t));
E = cell(N+1,numel(t));

hwb = waitbar(0,'Doing Matrix Multiplication for S');
for m = 0 : N
    for it = 1 : numel(t)
    Ecurrent = diag(1i.^(0:m));
    Gcurrent = diag(exp(1i.*(m:-2:-m)*t(it)));
    S0current = S0{m+1};
    Scurrent = real(Ecurrent*S0current'*Gcurrent*S0current*conj(Ecurrent));
    S{m+1,it} = Scurrent;
    G{m+1,it} = Gcurrent;
    E{m+1,it} = Ecurrent;
    end
    waitbar(m/N)
end
close(hwb)