function [ax_x,ax_y,ax_z] = fill_thi_axes(N,a,b,c)
%function [ax_x,ax_y,ax_z] = fill_thi_axes(N,a,b,c)

c0 = herm_coeff(0);
g = sqrt(a^2+b^2+c^2);
%analytic starting condition
thi_000 = c0^3*sqrt(2*pi)/g;


thi_00k = zeros(N+1,1); %z axis
thi_0k0 = zeros(N+1,1); %y axis
thi_k00 = zeros(N+1,1); %x axis

thi_00k(1) = thi_000; %starting value z
thi_0k0(1) = thi_000; %starting value y
thi_k00(1) = thi_000; %starting value x

i_even = 0:2:N;


for ii = 2 : numel(i_even)
    i_even_current = i_even(ii);
    i_even_prev = i_even(ii-1);
    k = i_even_current - 1;    
    %propagate z
    thi_00k(1+i_even_current) = z_flop_1d(k,thi_00k(1+i_even_prev),a,b,c); %we can use z flop as normal 
    %propate y
    thi_0k0(1+i_even_current) = z_flop_1d(k,thi_0k0(1+i_even_prev),a,c,b); %swap c for b   
    %propagate x 
    thi_k00(1+i_even_current) = z_flop_1d(k,thi_k00(1+i_even_prev),c,b,a); %swap c for a
       
    
end

ax_x = thi_k00;
ax_y = thi_0k0;
ax_z = thi_00k;

