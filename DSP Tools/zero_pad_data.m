function out = zero_pad_data(in,nzeros)
%function out = zero_pad_data(in,nzeros)

if numel(size(in))==2 && min(size(in))==1 %verify fbar is 1d
    data = in;
    data = [zeros(nzeros/2,1);data;zeros(nzeros/2,1)];
    out = data;
elseif numel(size(in))==2 && min(size(in))~=1 %verify fbar is 1d
    s1 = size(in,1); %assumes data is square
    s2 = size(in,2);
    if numel(nzeros) ==1
        nzeros1 = nzeros;
        nzeros2 = nzeros;
    elseif numel(nzeros)==2
        nzeros1 = nzeros(1);
        nzeros2 = nzeros(2);
    end  
    data = zeros(s1+nzeros1,s2+nzeros2);
    snew1 = size(data,1);
    snew2 = size(data,2);
    data(1+nzeros1/2:snew1-nzeros1/2,1+nzeros2/2:snew2-nzeros2/2) = in;
    out = data;
end