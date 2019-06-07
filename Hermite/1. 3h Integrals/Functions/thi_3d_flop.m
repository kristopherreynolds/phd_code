function [h_p1] = thi_3d_flop(m,n,k,a,b,c,h_m1)
%function h_p1 = thi_3d_flop(m,n,k,a,b,c,h_m1)
%
%
g = (a^2+b^2+c^2);

% M = diag([g/(2*a),g/(2*b),g/(2*c)]);
% N = [m*(2*a-g/a),2*b*n,2*c*k;2*a*m,n*(2*b-g/b),2*c*k;2*a*m,2*b*n,k*(2*c-g/c)];


D_H = [2*a*m*(2*a-g/a) 4*a*b*n 4*a*c*k; 4*a*b*m 2*b*n*(2*b-g/b) 4*b*c*k;...
    4*a*c*m 4*b*c*n 2*c*k*(2*c-g/c)]/g;

% D_H = M\N;

% cmp1 = herm_coeff(m+1);
% cmm1 = herm_coeff(m-1);
% cm = herm_coeff(m);
% cnp1 = herm_coeff(n+1);
% cnm1 = herm_coeff(n-1);
% cn = herm_coeff(n);
% ckp1 = herm_coeff(k+1);
% ckm1 = herm_coeff(k-1);
% ck = herm_coeff(k);
% 
% L_p1 = diag([cmp1*cn*ck cm*cnp1*ck cm*cn*ckp1]);
% L_m1_inv = diag(1./[cmm1*cn*ck cm*cnm1*ck cm*cn*ckm1]);
% 
% Q = L_p1*D_H*L_m1_inv;

% a1 = cmp1*cn*ck;
% a2 = cm*cnp1*ck;
% a3 = cm*cn*ckp1;
% 
% b1 = 1/(cmm1*cn*ck);
% b2 = 1/(cm*cnm1*ck);
% b3 = 1/(cm*cn*ckm1);

Q = zeros(3,3);
% % Q(1,1) = a1*b1*D_H(1,1);
% % Q(1,2) = a1*b2*D_H(1,2);
% % Q(1,3) = a1*b3*D_H(1,3);
% % Q(2,1) = a2*b1*D_H(2,1);
% % Q(2,2) = a2*b2*D_H(2,2);
% % Q(2,3) = a2*b3*D_H(2,3);
% % Q(3,1) = a3*b1*D_H(3,1);
% % Q(3,2) = a3*b2*D_H(3,2);
% % Q(3,3) = a3*b3*D_H(3,3);
% 
Q(1,1) = 0.5*sqrt(1/(m*(m+1)))*D_H(1,1);
Q(1,2) = 0.5*sqrt(1/(n*(m+1)))*D_H(1,2);
Q(1,3) = 0.5*sqrt(1/(k*(m+1)))*D_H(1,3);
Q(2,1) = 0.5*sqrt(1/(m*(n+1)))*D_H(2,1);
Q(2,2) = 0.5*sqrt(1/(n*(n+1)))*D_H(2,2);
Q(2,3) = 0.5*sqrt(1/(k*(n+1)))*D_H(2,3);
Q(3,1) = 0.5*sqrt(1/(m*(k+1)))*D_H(3,1);
Q(3,2) = 0.5*sqrt(1/(n*(k+1)))*D_H(3,2);
Q(3,3) = 0.5*sqrt(1/(k*(k+1)))*D_H(3,3);


h_p1 = Q*h_m1;

% %evaluate each new p1 value for the zero condition
% m0 = m;
% n0 = n;
% k0 = k;
% 
% %(m+1,n,k)
% m = m0+1;
% n = n0;
% k = k0;
% 
% sum_mn = (m+n)-k;
% sum_nk = (n+k)-m;
% sum_mk = (m+k)-n;
% evens_bool = ((m+n+k)/2 == round((m+n+k)/2));
% sum_bool1 = (sum_mn < 0);
% sum_bool2 = (sum_nk < 0);
% sum_bool3 = (sum_mk < 0);
% sum_bool = any([sum_bool1 sum_bool2 sum_bool3]);
% if  evens_bool && ~sum_bool
% 
% else    
% h_p1(1) = 0;    
% end
% 
% %(m,n+1,k)
% m = m0;
% n = n0+1;
% k = k0;
% 
% sum_mn = (m+n)-k;
% sum_nk = (n+k)-m;
% sum_mk = (m+k)-n;
% evens_bool = ((m+n+k)/2 == round((m+n+k)/2));
% sum_bool1 = (sum_mn < 0);
% sum_bool2 = (sum_nk < 0);
% sum_bool3 = (sum_mk < 0);
% sum_bool = any([sum_bool1 sum_bool2 sum_bool3]);
% if  evens_bool && ~sum_bool
% 
% else    
% h_p1(2) = 0;    
% end
% 
% %(m,n,k+1)
% m = m0;
% n = n0;
% k = k0+1;
% 
% sum_mn = (m+n)-k;
% sum_nk = (n+k)-m;
% sum_mk = (m+k)-n;
% evens_bool = ((m+n+k)/2 == round((m+n+k)/2));
% sum_bool1 = (sum_mn < 0);
% sum_bool2 = (sum_nk < 0);
% sum_bool3 = (sum_mk < 0);
% sum_bool = any([sum_bool1 sum_bool2 sum_bool3]);
% if  evens_bool && ~sum_bool
% 
% else    
% h_p1(3) = 0;    
% end
















