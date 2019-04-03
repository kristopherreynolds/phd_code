function h_p1 = xy_flop_2d(m,n,a,b,c,h_m1)
%function h_p1 = xy_flop2d(m,n,a,b,c,h_m1)
%
%
g = (a^2+b^2+c^2);

M = diag([g/(2*a) g/(2*b)]);
N = [m*(2*a-g/a) 2*b*n;2*a*m n*(2*b-g/b)];


%D_H = M\N;

D_H = [-(2*m*(-2*a^2+g))/g,4*a*b*n/g;(4*a*b*m)/g,-(2*(-2*n*b^2+g*n))/g];

cmp1 = herm_coeff(m+1);
cmm1 = herm_coeff(m-1);
cm = herm_coeff(m);
cnp1 = herm_coeff(n+1);
cnm1 = herm_coeff(n-1);
cn = herm_coeff(n);



L_p1 = diag([cmp1*cn cm*cnp1]);
L_m1 = diag([cmm1*cn cm*cnm1]);
L_m1inv = diag(1./(diag(L_m1)));

d11 = D_H(1,1);
d12 = D_H(1,2);
d21 = D_H(2,1);
d22 = D_H(2,2);

a1 = L_p1(1,1);
a2 = L_p1(2,2);
b1 = L_m1inv(1,1);
b2 = L_m1inv(2,2);

Q =  L_p1*D_H*L_m1inv;

% Q = [a1*b1*d11,a1*b2*d12;a2*b1*d21,a2*b2*d22];
% 
% Q = [cmp1*d11/cmm1,cmp1*cn*d12/(cm*cnm1);cm*cnp1*d21/(cmm1*cn),cnp1*d22/cnm1];
% 
Q(1,1) = d11*2^((m-1)/2-(m+1)/2)*sqrt(1/(m*(m+1)));
Q(2,2) = d22*2^((n-1)/2-(n+1)/2)*sqrt(1/(n*(n+1)));
Q(1,2) = 0.5*sqrt(1/(n*(m+1)))*d12;
Q(2,1) = 0.5*sqrt(1/(m*(n+1)))*d21;


h_p1 =Q*h_m1;

if cn < eps
    %keyboard
end

%replace nans
h_p1(isnan(h_p1)) = 0;


