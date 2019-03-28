%1d test script
clear
close all
clc




%create a square characteristic function
x = linspace(-5,5,2^8);
y = x;
step_width = 2;
I = zeros(numel(x),numel(y));
for ii = 1 : numel(x)
    for jj = 1 : numel(y)
        if x(ii) >=-step_width/2 && x(ii)<=step_width/2 && y(jj)>=-step_width/2 && y(jj)<=step_width/2
        I(ii,jj) = 1;
        end
    end
end
%zero pad
nzeros = [size(I,1)];
data = zero_pad_data(I,nzeros);

%load peppers image
I0 = imread('peppers.png');
I = rgb2gray(I0);
I = double(I);
nzeros = [512*2-384,512];  
data = zero_pad_data(I,nzeros);

clearvars -except data


P.N = 420;
P.range =28;
P.estimation_method = 'svd';
P.clims = [0 255];
P.cmap = 'gray';
P.special_limits = 0;
P.a = 1;
P.plot_error = 1;
hobj = hermite(P,data);
hobj.fit_series
hobj.plot

hobj
