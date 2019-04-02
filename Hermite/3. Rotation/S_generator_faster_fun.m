function S = S_generator_faster_fun(N,t)
S = cell(1,N+1);
hwb = waitbar(0,'Making S^m: Please Wait...');
for m = 0 : N
    Omega = zeros(m+1,m+1);
    for q = 0 : m
        for n = 0 : m
            Omega(q+1,n+1) = preSmatrix(q,n,m);
        end
    end
    S{m+1} = expm(Omega*t);
    waitbar(m/N)
end
close(hwb)
