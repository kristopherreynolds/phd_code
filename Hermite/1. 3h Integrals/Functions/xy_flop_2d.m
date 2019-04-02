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

h_p1 = L_p1*D_H/(L_m1)*h_m1;

%replace nans
h_p1(isnan(h_p1)) = 0;


