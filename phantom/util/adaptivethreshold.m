function bw=adaptivethreshold(I,ws,C)
bw = zeros( size(I), class(I));
[n, m] = size(I);

for i = 1:ws:n
    for j = 1:ws:m
        
        J = I(i:min(i+ws,n),j:min(j+ws,m));
        
        if max( J(:) > C) 
            level = graythresh( J );
            bw(i:min(i+ws,n),j:min(j+ws,m)) = J > level ; 
        end
    end
end