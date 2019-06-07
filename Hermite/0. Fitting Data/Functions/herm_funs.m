function out = herm_funs(N,x)
%function out = herm_funs(N,x)
% creates hermite functions from k = 0 to N for a vector x
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
    H = h0;
elseif N==1
    H = zeros(N+1,length(x));
    H(1,:) = h0;
    H(2,:) = h1;
elseif N==2
    H = zeros(N+1,length(x));
    H(1,:) = h0;
    H(2,:) = h1;
    H(3,:) = h2;
elseif N ==3
    H = zeros(N+1,length(x));
    H(1,:) = h0;
    H(2,:) = h1;
    H(3,:) = h2;
    H(4,:) = h3;
else
    H = zeros(N+1,length(x));
    H(1,:) = h0;
    H(2,:) = h1;
    H(3,:) = h2;
    H(4,:) = h3;
    for k = 4:N
        n = k-1;
        hn = H(k,:);
        hnm1 = H(k-1,:);
        h4 = ((x*sqrt(2)).*hn - sqrt(n)*hnm1)/sqrt(n+1);
        H(k+1,:) = h4;
    end
end
out = H';

end

