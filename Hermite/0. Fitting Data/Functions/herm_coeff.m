function out = herm_coeff(n)
%function cn = herm_coeff(n)
out = 1./(2.^(n./2).*sqrt(factorial(n).*sqrt(pi)));
end
