function xs = completePathImageSequence( xs, x0, x1, flipSequence )

k = size(xs, 3);

if flipSequence
    xt = x0;
    x0 = x1;
    x1 = xt;
end



j = 0;
for i = 1 : k
    xt = xs(:,:,i);
    if max( xt(:)) < 0.01 
        j = j + 1;
    end
end


if j > 0
    xs(:,:,2:end) = xs(:,:,1:end-1);
    xs(:,:,1) = x0;
    j = j - 1; 
end


if j > 0
    xs(:,:,end) = x1;
    j = j - 1;
end

if j > 0
    for i = 1:j
        
        xs(:,:,end-i) = (i/(j+1)) * xs(:,:,end) + ((j-i+1)/(j+1)) * xs(:,:,end-j-1);
        
        
    end
end


if flipSequence
    xs(:,:,:) = xs(:,:,end:-1:1);
end



end