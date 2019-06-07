%test z flop 1d

clear
close all
clc

N = 10; %band limit

a = 2;
b = a;
c = sqrt(a^2+b^2);
g = sqrt(a^2+b^2+c^2);

c0 = herm_coeff(0);

%analytic starting conditions
thi_000 = c0^3*sqrt(2*pi)/g;

[thi_k00,thi_0k0,thi_00k] = fill_thi_axes(N,a,b,c);

thi_k00anl = zeros(N+1,1);
thi_0k0anl = zeros(N+1,1);
thi_00kanl = zeros(N+1,1);
xpts = linspace(-15,15,2000);

est_type = 'anl';
for k = 0 : N
    
    if strcmp(est_type,'num')
    thi_k00anl(k+1) = thi_num(k,0,0,a,b,c,xpts);
    thi_0k0anl(k+1) = thi_num(0,k,0,a,b,c,xpts);
    thi_00kanl(k+1) = thi_num(0,0,k,a,b,c,xpts);
    elseif strcmp(est_type,'anl')
    thi_k00anl(k+1) = thi_anl2(k,0,0,a,b,c,2);
    thi_0k0anl(k+1) = thi_anl2(0,k,0,a,b,c,2);
    thi_00kanl(k+1) = thi_anl2(0,0,k,a,b,c,2);
    end
    
end

err_vec_z = abs(thi_00k-thi_00kanl);
err_vec_y = abs(thi_0k0-thi_0k0anl);
err_vec_x = abs(thi_k00-thi_k00anl);

figure
subplot(3,1,1)
plot(0:N,thi_k00,0:N,thi_k00anl,'m--','linewidth',2)
title('rr')
grid on
legend('RR','True')
subplot(3,1,2)
plot(0:N,thi_0k0,0:N,thi_0k0anl,'m--','linewidth',2)
grid on
legend('RR','True')
subplot(3,1,3)
plot(0:N,thi_00k,0:N,thi_00kanl,'m--','linewidth',2)
grid on
xlabel('k')
legend('RR','True')







