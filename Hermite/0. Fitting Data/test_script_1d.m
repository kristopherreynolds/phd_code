%1d test script
clear
close all
clc




%create a step function
x = linspace(-5,5,2^8);
y = x*0;

step_width =2;
i1 = x>=-step_width/2;
i2 = x<=step_width/2;
i12 = find(i1==1&i2==1);
y(i12) = 1;

data = y';clearvars -except data
nzeros = 2^8;
%zero pad data
data = zero_pad_data(data,nzeros);

P.N = 420;
P.range =28;
P.estimation_method = 'svd';
P.a = 1;
P.plot_error = 1;
hobj = hermite(P,data);
hobj.fit_series
hobj.plot

hobj
