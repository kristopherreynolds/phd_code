function [Sxy,Syz,Sxz] = fill_thi_planes(N,a,b,c)
%function [Sxy,Syz,Sxz] = fill_thi_planes(N,a,b,c)


Sxy = zeros(N+1,N+1);
Syz = zeros(N+1,N+1);
Sxz = zeros(N+1,N+1);


[ax_x,ax_y,ax_z] = fill_thi_axes(N,a,b,c); %uses 1d recurrence relation

%keystone points
c1 = herm_coeff(1);
c0 = herm_coeff(0);

thi_110 = c1^2*c0*4*sqrt(2*pi)*a*b/((a^2+b^2+c^2)^(3/2));
thi_011 = c1^2*c0*4*sqrt(2*pi)*c*b/((a^2+b^2+c^2)^(3/2)); 
thi_101 = c1^2*c0*4*sqrt(2*pi)*a*c/((a^2+b^2+c^2)^(3/2)); 

%initialize xy plane
Sxy(:,1) = ax_x;
Sxy(1,:) = ax_y;
Sxy(2,2) = thi_110;

%initialize yz plane
Syz(:,1) = ax_y;
Syz(1,:) = ax_z;
Syz(2,2) = thi_011;

%initialize xz plane
Sxz(:,1) = ax_x;
Sxz(1,:) = ax_z;
Sxz(2,2) = thi_101;

M = 2e3; %large arbitary number of iterations. for loop has a break condition that will proc first

%initialize ref points
ref_pts = [1 1];
for ipass = 1 : M
    
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
        disp('done')
        break
    end
    
    %loop over ref points
    for ipts = 1 : size(ref_pts,1)
        ref_point_current = ref_pts(ipts,:);
        
        m = ref_point_current(:,1);
        n = ref_point_current(:,2);        
        
        %apply RR for xy plane
        thiold1_xy = Sxy(m-1+1,n+1);
        thiold2_xy = Sxy(m+1,n-1+1);
        thi_m1_xy = [thiold1_xy;thiold2_xy];
        
        if sum(thi_m1_xy(:))~=0
        thi_p1_xy = xy_flop_2d(m,n,a,b,c,thi_m1_xy);
        else
        thi_p1_xy = thi_m1_xy;
        end
        
        pairs = [m+1,n;m,n+1];
        odds = sum(pairs,2);
        odds_bool = (round(odds/2)~=odds/2);
        thi_p1_xy(odds_bool==1) = 0;
        
        Sxy(m+1+1,n+1) = thi_p1_xy(1);
        Sxy(m+1,n+1+1) = thi_p1_xy(2);
        
        %apply RR for yz plane (swap c for a)
        thiold1_yz = Syz(m-1+1,n+1);
        thiold2_yz = Syz(m+1,n-1+1);
        thi_m1_yz = [thiold1_yz;thiold2_yz];
        
        if sum(thi_m1_yz(:))~=0
            thi_p1_yz = xy_flop_2d(m,n,b,c,a,thi_m1_yz);
        else
            thi_p1_yz = thi_m1_yz;
        end
        
        pairs = [m+1,n;m,n+1];
        odds = sum(pairs,2);
        odds_bool = (round(odds/2)~=odds/2);
        thi_p1_yz(odds_bool==1) = 0;
        
        Syz(m+1+1,n+1) = thi_p1_yz(1);
        Syz(m+1,n+1+1) = thi_p1_yz(2);
        
        %apply RR for xz plane (swap c for b)
        thiold1_xz = Sxz(m-1+1,n+1);
        thiold2_xz = Sxz(m+1,n-1+1);
        thi_m1_xz = [thiold1_xz;thiold2_xz];
        
        if sum(thi_m1_xz(:))~=0
            thi_p1_xz = xy_flop_2d(m,n,a,c,b,thi_m1_xz);
        else
            thi_p1_xz = thi_m1_xz;
        end
        
        pairs = [m+1,n;m,n+1];
        odds = sum(pairs,2);
        odds_bool = (round(odds/2)~=odds/2);
        thi_p1_xz(odds_bool==1) = 0;
        
        Sxz(m+1+1,n+1) = thi_p1_xz(1);
        Sxz(m+1,n+1+1) = thi_p1_xz(2);
        
        
    end
    
    %after looping over all new points, update ref_points
    ref_pts = newref_pts;
    
end

%Crop
Sxy = Sxy(1:N+1,1:N+1);
Syz = Syz(1:N+1,1:N+1);
Sxz = Sxz(1:N+1,1:N+1);









