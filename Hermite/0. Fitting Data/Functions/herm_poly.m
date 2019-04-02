function out = herm_poly(n,x)
%out = herm_poly(n,x)
%version 2
%first, choose the appropriate Hermite Polynomial given the order 'n'
if n == 0
    H = 1*ones(size(x));
elseif n == 1
    H = 2.*x;
elseif n == 2
    H = 4.*x.^2 - 2*ones(size(x));
elseif n == 3
    H = 8.*x.^3 - 12.*x;
elseif n > 3 %use recursion to get higher order polynomials
    diff = n-3; %the number of times to apply recursion
    Hn  = 8.*x.^3 - 12.*x; 
    Hnm1 = 4.*x.^2 - 2*ones(size(x));  
    n0 = 3;
    %apply recurrence formula
    for ii = 1 : diff
    Hnp1 = 2*x.*Hn - 2*(n0+(ii-1))*Hnm1;
    %reset the previous iterations
    Hnm1 = Hn;
    Hn = Hnp1;
    end
    H = Hnp1;   
end
out = H;