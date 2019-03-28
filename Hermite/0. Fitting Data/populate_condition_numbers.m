%Populate Condition Numbers

npts = 2*512;
Nvec = 300:500;
Rvec = linspace(25,30,25);
C = condition_number(Rvec,Nvec,npts);
figure;
imagesc(Nvec,Rvec,log(abs(C)))