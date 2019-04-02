function h_p1 = thi_3d_flop(m,n,k,a,b,c,h_m1)
%function h_p1 = thi_3d_flop(m,n,k,a,b,c,h_m1)
%
%
g = (a^2+b^2+c^2);

% M = diag([g/(2*a) g/(2*b) g/(2*c)]);
% N = [m*(2*a-g/a) 2*b*n 2*c*k;2*a*m n*(2*b-g/b) 2*c*k;2*a*m 2*b*n k*(2*c-g/c)];

D_H = [2*a*m*(2*a-g/a) 4*a*b*n 4*a*c*k; 4*a*b*m 2*b*n*(2*b-g/b) 4*b*c*k;...
    4*a*c*m 4*b*c*n 2*c*k*(2*c-g/c)]/g;

cmp1 = herm_coeff(m+1);
cmm1 = herm_coeff(m-1);
cm = herm_coeff(m);
cnp1 = herm_coeff(n+1);
cnm1 = herm_coeff(n-1);
cn = herm_coeff(n);
ckp1 = herm_coeff(k+1);
ckm1 = herm_coeff(k-1);
ck = herm_coeff(k);

L_p1 = diag([cmp1*cn*ck cm*cnp1*ck cm*cn*ckp1]);
L_m1 = diag([cmm1*cn*ck cm*cnm1*ck cm*cn*ckm1]);

h_p1 = L_p1*D_H*inv(L_m1)*h_m1;
