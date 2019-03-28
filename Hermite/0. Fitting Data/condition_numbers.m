function out = condition_numbers(Rvec,Nvec,npts)

C = zeros(numel(Rvec),numel(Nvec));
inds = 1 : npts;
hwb = waitbar(0,'wait');
for ii = 1 : numel(Rvec)
    for jj = 1 : numel(Nvec)
        Rcurrent = Rvec(ii);
        x = -Rcurrent+2*Rcurrent*(inds-1)/(npts-1);
        H = herm_funs(Nvec(jj),x);
        [~,S,~] = svd(H);
        s = diag(S);
        s = s(1:Nvec(jj)+1);
        C(ii,jj) = max(s)/min(s);
        
    end
    waitbar(ii/numel(Rvec))
end
close(hwb)

out = C;