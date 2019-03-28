function S = S_generator_fun(mmax,theta)
%function S = S_generator_fun(mmax,theta)

%S Matrix Generator

%% Initialize S by filling in S^0, S^1, and S^2

S{0+1} = 1*size(theta);
S{1+1} = [cos(theta) sin(theta);-sin(theta) cos(theta)];
S{2+1} = [cos(theta).^2 sqrt(2)*cos(theta).*sin(theta) sin(theta).^2;-sqrt(2)*cos(theta).*sin(theta) cos(theta).^2-sin(theta).^2 sqrt(2)*cos(theta).*sin(theta);...
    sin(theta).^2 -sqrt(2).*cos(theta).*sin(theta) cos(theta).^2];

%% Generate S^m 

hwb = waitbar(0,'Making S^m: Please Wait...');
clc
for m = 2 : mmax-1
    S = propagate_S(S,m,theta);
    waitbar(m/(mmax-1))
end
close(hwb)
 
 S{1} = 1;