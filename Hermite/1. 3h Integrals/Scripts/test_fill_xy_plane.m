%test fill xy plane

%PLANAR RECURRENCE RELATION ILLUSTRATION FOR 3h INTEGRALS
%Author: KLR
%Description: This code will illustrate and produce the results for generating
% a N X N integter block square of results for corresponding 3h integrals
% based off of their Hermite function order. The coordinates for a 3h
% integral are m,n, and k (k=0 for planar). m is x-axis and n is y-axis

%Date Last Modified: 3/3/2018


%% Clear Workspace
clear
close all
%clc

%% Define Initial Conditions

fsz = 14; %fontsize

N = 200; %the limits of NxN square
M = 2e3; %large arbitary number of iterations

a = 1;
b = 1;
c = 1;

%x-axis points (m-axis)
x = 0:N;
y = 0:N*0;
[X,Y] = meshgrid(x,y);
xv1 = X(:);
yv1 = Y(:); clear x y X Y

%y-axis points (n-axis)
x = 0:N*0;
y = 0:N;
[X,Y] = meshgrid(x,y);
xv2 = X(:);
yv2 = Y(:); clear x y X Y

%add keystone point
keystone = [1 1];
xv = [xv1;xv2;keystone(1)];
yv = [yv1;yv2;keystone(2)];

%store initial conditions in a set variable called 'pt_set' that
%will constantly change in size
pt_set = [xv yv];


%% Set up Figure with starting 3h integrals

%now, lets populate rr version of S, filling it with the starting
%conditions
Srr = zeros(N+3,N+3);
Sanl = zeros(N+3,N+3);

[ax_x,ax_y,~] = fill_thi_axes(N+2,a,b,c); %uses 1d recurrence relation
Srr(:,1) = ax_x;
Srr(1,:) = ax_y;
Sanl(:,1) = ax_x;
Sanl(1,:) = ax_y;
c1 = herm_coeff(1);
c0 = herm_coeff(0);
thi_110 = c1^2*c0*4*sqrt(2*pi)*a*b/((a^2+b^2+c^2)^(3/2));
Srr(2,2) = thi_110;
Sanl(2,2) = thi_110;


ref_pts = keystone; %keystone is starting reference point

%% Now for the full RR Propagation
tic

%perform passes
for ipass = 1 : M
    
    %create new reference points by selecting current ones and adding +1 to either m
    %or n vector
    %mvec = ref_pts(:,1);
    %nvec = ref_pts(:,2);
    %newref_pts = unique([mvec+1 nvec;mvec nvec+1],'rows');
    %c2 = newref_pts(:,2);
    %[~,isort] = sort(c2); clear c3
    %newref_pts = newref_pts(isort,:); clear isort
    newref_pts = [(1:ipass+1)',(ipass+1:-1:1)'];

    %now, loop over reference points, and eliminate any that have
    %(m,n,k) values greater than N
    maxes = zeros(size(newref_pts,1),1);
    
    for ipts = 1 : size(newref_pts,1)
        pt_current = newref_pts(ipts,:);
        [ymax,~] = max(pt_current);
        maxes(ipts) = ymax;
    end
    newref_pts = newref_pts(maxes<=N,:);
    if isempty(newref_pts)
        npasses = ipass-1;
        break
    end
    
    %title(['Pass #',num2str(ipass)])
    %disp(['Pass # ',num2str(ipass)])
    
    %loop over ref points
    for ipts = 1 : size(ref_pts,1)
        ref_point_current = ref_pts(ipts,:);
        
        m = ref_point_current(:,1);
        n = ref_point_current(:,2);
        old_pts = [m-1 n;m n-1];
        new_pts = [m+1 n;m n+1];
        
        
        %apply RR
        thiold1 = Srr(m-1+1,n+1);
        thiold2 = Srr(m+1,n-1+1);
        thi_m1 = [thiold1;thiold2];
        if sum(thi_m1(:))~=0
        thi_p1 = xy_flop_2d(m,n,a,b,c,thi_m1);
        else
        thi_p1 = thi_m1;
        end

        Srr(m+1+1,n+1) = thi_p1(1);
        Srr(m+1,n+1+1) = thi_p1(2);
        
    end
    
    %after looped over all new points, update ref_points
    ref_pts = newref_pts;
    
end

timeout = toc;

%crop
Srr = Srr(1:N+1,1:N+1);

Snum = Srr*0;
tic
x = linspace(-15,15,1000);
for ii = 0 : N
    for jj = 0 : N
        Snum(ii+1,jj+1) = thi_num(ii,jj,0,a,b,c,x);
    end
end
timeout_num = toc;



subplot(1,2,1)
imagesc(Srr)
axis xy
colorbar
caxis([-0.2 0.2])
title('RR')
subplot(1,2,2)
imagesc(Snum)
axis xy
colorbar
caxis([-0.2 0.2])
title('Num')

l2_err = norm(Srr-Snum);


% figure
% imagesc(0:N+2,0:N+2,abs(Srr-Sanl)')
% axis xy
% cb = colorbar;
% ylabel(cb,'Absolute Difference')
% xlabel('m')
% ylabel('n')
% title('3h Integral Error from Recurrence Relations')
% set(gca,'xlim',[0 7],'ylim',[0 7],'fontsize',fsz)
% hold on
% plot([N N],[-0.5 N],'k')
% plot([-0.5 N],[N N],'k')
% caxis([0 5]*1e-16)


% Snum = Srr*0;
% xpts = linspace(-40,40,50000);
% for m = 0 : N
%     for n = 0 : N
%         Snum(m+1,n+1) = thi_num(m,n,0,a,b,c,xpts);
%     end
% end
% 
% figure
% imagesc(0:N+2,0:N+2,abs(Snum-Sanl)')
% axis xy
% cb = colorbar;
% ylabel(cb,'Absolute Difference')
% xlabel('m')
% ylabel('n')
% title('3h Integral Error from Numeric Integration')
% set(gca,'xlim',[0 7],'ylim',[0 7],'fontsize',fsz)
% hold on
% plot([N N],[-0.5 N],'k')
% plot([-0.5 N],[N N],'k')
% caxis([0 5]*1e-16)






