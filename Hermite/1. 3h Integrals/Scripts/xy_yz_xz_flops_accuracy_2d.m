clear
close all
clc

N = 150;

a = 1;
b = 1;
c = sqrt(2);


% Sxy = fill_xy_plane(N,a,b,c);
% Syz = fill_xy_plane(N,b,c,a);
% Sxz = fill_xy_plane(N,a,c,b);

[Sxy,Syz,Sxz] = fill_thi_planes(N,a,b,c);

Sxynum = Sxy*0;
Syznum = Syz*0;
Sxznum = Sxz*0;
x = linspace(-15,15,1000);
for ii = 0 : N
    for jj = 0 : N
        Sxynum(ii+1,jj+1) = thi_num(ii,jj,0,a,b,c,x);
        Syznum(ii+1,jj+1) = thi_num(0,ii,jj,a,b,c,x);
        Sxznum(ii+1,jj+1) = thi_num(ii,0,jj,a,b,c,x);
    end
end

fsz = 14;

%XY
figure
subplot(1,2,1)
imagesc(Sxy)
xlabel('X')
ylabel('Y')
axis xy
colorbar
title('RR')
set(gca,'fontsize',fsz)
subplot(1,2,2)
imagesc(abs(Sxy-Sxynum))
axis xy
colorbar
title('|RR-Num|')
set(gca,'fontsize',fsz)
xlabel('X')
ylabel('Y')


%Yz
figure
subplot(1,2,1)
imagesc(Syznum)
xlabel('Y')
ylabel('Z')
axis xy
colorbar
title('RR')
set(gca,'fontsize',fsz)
subplot(1,2,2)
imagesc(abs(Syz-Syznum))
axis xy
colorbar
title('|RR-Num|')
set(gca,'fontsize',fsz)
xlabel('Y')
ylabel('Z')


%XZ
figure
subplot(1,2,1)
imagesc(Sxznum)
xlabel('X')
ylabel('Z')
axis xy
colorbar
title('RR')
set(gca,'fontsize',fsz)
subplot(1,2,2)
imagesc(abs(Sxz-Sxznum))
axis xy
colorbar
title('|RR-Num|')
set(gca,'fontsize',fsz)
xlabel('X')
ylabel('Z')