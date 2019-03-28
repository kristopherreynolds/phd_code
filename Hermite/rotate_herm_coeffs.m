function fhat_rot = rotate_herm_coeffs(fhat,S,N)
%function fhat_rot = rotate_herm_coeffs(fhat,S,N)
%
% fhat is a matrix of size (N+1)x(N+1)
% fmt = Sm*fm, where 
% fm = ..

fhat = squeeze(fhat);
%compute rotated coefficients
fhat_rot = zeros(size(fhat));
fhat_rot = fhat;

for m = 0 : N   
    
%     if m < N
%         break
%     end

    %z rotation
    i2 = m:-1:0;
    i1 = 0:1:m;
    
    %x rotation
    
    fm = diag(fhat(1+i1,1+i2));  
    Sm = S{m+1}(:,:);
    fmt = Sm*fm;    
    idx = sub2ind(size(fhat_rot),1+i1,1+i2);
    fhat_rot(idx) = fliplr(fmt);
    
    %[fm fliplr(fmt)]
%     dfm = fm-fliplr(fmt);
%     plot(dfm)
%     keyboard
    
end





