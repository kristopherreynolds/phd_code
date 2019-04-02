%test z flop 1d

clear
close all
clc

N = 550; %band limit

a = 1;
b = 2;
c = 3;
g = sqrt(a^2+b^2+c^2);

c0 = herm_coeff(0);

%analytic starting conditions
thi_000 = c0^3*sqrt(2*pi)/g;

%because of the 3 results below, all odd flops are zero too since they are derived from this
thi_001 = 0; 
thi_010 = 0;
thi_100 = 0;


i_even = 0:2:N;

method = 2; %symbolic integration method (select 1 or 2)

int_type = 'num';

thi_00k = zeros(N+1,1); %z
thi_0k0 = zeros(N+1,1); %y
thi_k00 = zeros(N+1,1); %x
thi_00kanl = zeros(N+1,1); %zanl
thi_0k0anl = zeros(N+1,1); %yanl
thi_k00anl = zeros(N+1,1); %xanl

thi_00k(1) = thi_000; %starting value z
thi_0k0(1) = thi_000; %starting value y
thi_k00(1) = thi_000; %starting value x

if strcmp(int_type,'sym')
q = thi_anl(0,0,0,a,b,c,method);
thi_00kanl(1) = q.thi; clear q
q = thi_anl(0,0,0,a,c,b,method);
thi_0k0anl(1) = q.thi; clear q
q = thi_anl(0,0,0,c,b,a,method);
thi_k00anl(1) = q.thi; clear q
elseif strcmp(int_type,'num')
x = linspace(-15,15,1000);
dx = mean(diff(x));
thi_000num = thi_num(0,0,0,a,b,c,x);
thi_00kanl(1) = thi_000num;
thi_0k0anl(1) = thi_000num;
thi_k00anl(1) = thi_000num;
end

for ii = 2 : numel(i_even)
    i_even_current = i_even(ii);
    i_even_prev = i_even(ii-1);
    k = i_even_current - 1;
    
    %z
    thi_00k(1+i_even_current) = z_flop_1d(k,thi_00k(1+i_even_prev),a,b,c);
    if strcmp(int_type,'sym')
    q = thi_anl(0,0,i_even_current,a,b,c,method);
    thi_00kanl(1+i_even_current) = q.thi; clear q
    elseif strcmp(int_type,'num')
    thi_00kanl(1+i_even_current)= thi_num(0,0,i_even_current,a,b,c,x);    
    end
    %y
    thi_0k0(1+i_even_current) = z_flop_1d(k,thi_0k0(1+i_even_prev),a,c,b);
    if strcmp(int_type,'sym')
    q = thi_anl(0,i_even_current,0,a,b,c,method);
    thi_0k0anl(1+i_even_current) = q.thi; clear q
    elseif strcmp(int_type,'num')
    thi_0k0anl(1+i_even_current)= thi_num(0,i_even_current,0,a,b,c,x);    
    end
    %x
    thi_k00(1+i_even_current) = z_flop_1d(k,thi_k00(1+i_even_prev),c,b,a);
    if strcmp(int_type,'sym')
    q = thi_anl(i_even_current,0,0,a,b,c,method);
    thi_k00anl(1+i_even_current) = q.thi; clear q
    elseif strcmp(int_type,'num')
    thi_k00anl(1+i_even_current) = thi_num(i_even_current,0,0,a,b,c,x);    
    end
end

err_vec_z = abs(thi_00k-thi_00kanl);
err_vec_y = abs(thi_0k0-thi_0k0anl);
err_vec_x = abs(thi_k00-thi_k00anl);

figure
ax1 = subplot(2,1,1);
plot(0:N,thi_k00,0:N,thi_0k0,0:N,thi_00k,'linewidth',2)
legend('x','y','z')
grid on
xlabel('N')
ax2 = subplot(2,1,2);
plot(0:N,err_vec_x,0:N,err_vec_y,0:N,err_vec_z,'linewidth',1)
grid on
xlabel('N')
ylabel('Absolute Error')
legend('dx','dy','dz')
linkaxes([ax1 ax2],'x')

