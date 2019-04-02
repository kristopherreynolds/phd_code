%fill axes test script
clear
close all
clc

N = 540;
a = 1;
b = 2;
c = 3;

[ax_x,ax_y,ax_z] = fill_thi_axes(N,a,b,c);

figure
plot(0:N,ax_x,0:N,ax_y,0:N,ax_z);
legend('x','y','z')
xlabel('N')

