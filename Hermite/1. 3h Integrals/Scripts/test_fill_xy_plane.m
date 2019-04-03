


%% Clear Workspace
clear
close all
clc

%% Define Initial Conditions
clims = [-0.1 0.1];
fsz = 14; %fontsize

N = 150; %the limits of NxN square
M = 2e3; %large arbitary number of iterations

a = sqrt(2);
b = 1;
c = 1;


%% Fill out X & Y Axes

%now, lets populate rr version of S, filling it with the starting
%conditions
Srr = zeros(N+1,N+1);

[ax_x,ax_y,~] = fill_thi_axes(N,a,b,c); %uses 1d recurrence relation
Srr(:,1) = ax_x;
Srr(1,:) = ax_y;
c1 = herm_coeff(1);
c0 = herm_coeff(0);
thi_110 = c1^2*c0*4*sqrt(2*pi)*a*b/((a^2+b^2+c^2)^(3/2));
Srr(2,2) = thi_110;

% %add keystone point
keystone = [1 1];

ref_pts = keystone; %keystone is starting reference point

%% Fill XY Plane
% tic
% 
% %perform passes
% for ipass = 1 : M
%     
%     newref_pts = [(1:ipass+1)',(ipass+1:-1:1)'];
% 
%     %now, loop over reference points, and eliminate any that have
%     %(m,n,k) values greater than N
%     maxes = zeros(size(newref_pts,1),1);
%     
%     for ipts = 1 : size(newref_pts,1)
%         pt_current = newref_pts(ipts,:);
%         [ymax,~] = max(pt_current);
%         maxes(ipts) = ymax;
%     end
%     newref_pts = newref_pts(maxes<=N,:);
%     if isempty(newref_pts)
%         disp('done')
%         break
%     end
%     
%     %loop over ref points
%     for ipts = 1 : size(ref_pts,1)
%         ref_point_current = ref_pts(ipts,:);
%         
%         m = ref_point_current(:,1);
%         n = ref_point_current(:,2);        
%         
%         %apply RR
%         thiold1 = Srr(m-1+1,n+1);
%         thiold2 = Srr(m+1,n-1+1);
%         thi_m1 = [thiold1;thiold2];
%         if sum(thi_m1(:))~=0
%         thi_p1 = xy_flop_2d(m,n,a,b,c,thi_m1);
%         else
%         thi_p1 = thi_m1;
%         end
% 
%         Srr(m+1+1,n+1) = thi_p1(1);
%         Srr(m+1,n+1+1) = thi_p1(2);
%         
%     end
%     
%     %after looped over all new points, update ref_points
%     ref_pts = newref_pts;
%     
% end
% 
% timeout = toc;
% 
% %crop
% Srr = Srr(1:N+1,1:N+1);
% 
% 
% figure
% imagesc(0:N,0:N,db(Srr))
% caxis([-70 10])
% axis xy

Srr = fill_xy_plane(N,a,b,c);


%numeric

Snum = Srr*0;
tic
x = linspace(-15,15,1000);
for ii = 0 : N
    for jj = 0 : N
        Snum(ii+1,jj+1) = thi_num(ii,jj,0,a,b,c,x);
    end
end
timeout_num = toc;


figure
subplot(1,2,1)
imagesc(Snum)
xlabel('M')
ylabel('N')
axis xy
colorbar
caxis(clims)
title('RR')
set(gca,'fontsize',fsz)
subplot(1,2,2)
imagesc(abs(Srr-Snum))
axis xy
colorbar
title('|RR-Num|')
set(gca,'fontsize',fsz)
xlabel('M')
ylabel('N')

l2_err = norm(Srr-Snum);









