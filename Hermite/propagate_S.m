function S = propagate_S(S,m,theta)

for n = 0 :m
    for q = 0 : m
        
        out1 = S_rr_eq1(S,n,q,m,theta);
        out2 = S_rr_eq2(S,n,q,m,theta);
        
        if isstruct(out1) == 1
            %add values to make S3
            %first value
            ii = 1;
            iq_current = out1.indices(ii,1);
            in_current = out1.indices(ii,2);
            im_current = out1.indices(ii,3);
            S{im_current}(iq_current,in_current) = out1.svalues(ii);
            %second value
            ii = 2;
            iq_current = out1.indices(ii,1);
            in_current = out1.indices(ii,2);
            im_current = out1.indices(ii,3);
            S{im_current}(iq_current,in_current) = out1.svalues(ii);
        end
        
        if isstruct(out2) == 1
            %add values to make S3
            %first value
            ii = 1;
            iq_current = out2.indices(ii,1);
            in_current = out2.indices(ii,2);
            im_current = out2.indices(ii,3);
            S{im_current}(iq_current,in_current) = out2.svalues(ii);
            %second value
            ii = 2;
            iq_current = out2.indices(ii,1);
            in_current = out2.indices(ii,2);
            im_current = out2.indices(ii,3);
            S{im_current}(iq_current,in_current) = out2.svalues(ii);
        end
        
    end
end

%fill in corner cases
S{im_current}(0+1,0+1) = cos(theta)^(m+1);
S{im_current}(im_current,im_current) = cos(theta)^(m+1);
S{im_current}(0+1,im_current) = sin(theta)^(m+1);
S{im_current}(im_current,0+1) = (-sin(theta))^(m+1);
