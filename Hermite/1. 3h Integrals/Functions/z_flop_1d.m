function thi_k_plus_1 = z_flop_1d(k,thi_k_minus_1,a,b,c)
%function thi_k_plus_1 = z_flop_1d(k,thi_k_minus_1,a,b,c)

g = (a^2+b^2+c^2);

%ckp1 = herm_coeff(k+1);
%ckm1 = herm_coeff(k-1);

%c_ratio = ckp1/ckm1;
%more stable c_ratio
c_ratio = (1/2)*sqrt(1/(k*(k+1)));


thi_k_plus_1  = c_ratio*(2*k*c)/g*(2*c-g/c)*thi_k_minus_1;

% if k > 170
%     %keyboard
% end

if isnan(thi_k_plus_1)
    thi_k_plus_1 = 0;
end


