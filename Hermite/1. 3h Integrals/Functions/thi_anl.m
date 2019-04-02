function out = thi_anl(m,n,k,a,b,c,method)
%function out = thi_anl(m,n,k,a,b,c,method)
%out.tHi = int[Hm(ax)Hn(bx)Hk(cx)exp(-x^2/2)] from x = -inf to inf
%out.thi = cmcnck*int[Hm(ax)Hn(bx)Hk(cx)exp(-x^2/2)] from x = -inf to inf

syms x real


Hm = herm_poly(m,a*x);
Hn = herm_poly(n,b*x);
Hk = herm_poly(k,c*x);

g = (a^2+b^2+c^2);

out.tHi = eval(int(Hm*Hn*Hk*exp(-g*x^2/2),x,-inf,inf));

cm = herm_coeff(m);
cn = herm_coeff(n);
ck = herm_coeff(k);

if method ==1
out.thi = cm*cn*ck*out.tHi;
elseif method == 2
hm = herm_fun(m,a*x);
hn = herm_fun(n,b*x);
hk = herm_fun(k,c*x);
out.thi = eval(int(hm*hn*hk,x,-inf,inf));
end




