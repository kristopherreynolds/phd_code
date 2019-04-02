function h = herm_fun(N,x)
%h = herm_fun(n,x)

% H = herm_poly(n,x);
% cn = herm_coeff(n);
% h = cn*H.*exp(-x.^2/2);

%first, lets size up x. if it is a column vector, make it a row vector
s = size(x);
s1 = s(1);
s2 = s(2);
if any(s==1)
    if s1>s2
        x = x';
    end
end

h0 = herm_coeff(0).*exp(-x.^2/2);
h1 = herm_coeff(1).*(2*x).*exp(-x.^2/2);
h2 = herm_coeff(2).*(4.*x.^2 - 2*ones(size(x))).*exp(-x.^2/2);
h3 = herm_coeff(3).*(8.*x.^3 - 12.*x).*exp(-x.^2/2);

if N == 0
    h = h0;
elseif N==1
    h = h1;
elseif N==2
    h = h2;
elseif N ==3
    h = h3;
else
    hn = h3;
    hnm1 = h2;
    for k = 4:N
        n = k-1;
        hnp1 = ((x*sqrt(2)).*hn - sqrt(n)*hnm1)/sqrt(n+1);
        hnm1 = hn;
        hn = hnp1;
    end
    h = hnp1;
end


