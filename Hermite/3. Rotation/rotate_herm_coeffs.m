function fhat_rot = rotate_herm_coeffs(fhat,S,N)
%function fhat_rot = rotate_herm_coeffs(fhat,S,N)
%
% fhat is a matrix of size (N+1)x(N+1)
% fmt = Sm*fm, where 
% fm = ..

fhat = squeeze(fhat);
%compute rotated coefficients
fhat_rot = zeros(size(fhat));



for m = 0 : N   
    
    
    i1 = 0:1:m;
    i2 = m:-1:0;
    
    idx = sub2ind(size(fhat),1+i1,1+i2);
    fm = fhat(idx)';
     
    Sm = S{m+1};
    fmt = Sm*fm;    
    idx = sub2ind(size(fhat_rot),1+i1,1+i2);
    fhat_rot(idx) = fmt;
       
end





