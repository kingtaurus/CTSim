function J = minimizeLocalStd( I, C, n, map )

m = (n - 1)/2;

weight = zeros( size(I) );
J = zeros( size(I) );

for i = n : size(I, 1) - n
    for j = n : size( I, 2) - n

        X = I( i-m:i+m, j-m:j+m );
        Y = C( i-m:i+m, j-m:j+m );
        
        if  map(i,j) == 0 || max(abs(Y(:))) < 1e-4 
            continue;
        end
        
        
        
        X = X - mean(X(:));
        Z = Y - mean(Y(:));
        
        w = sum( X(:).* Z(:) ) / sqrt(sum( Z(:).*Z(:) )) ;
        
        J( i-m:i+m, j-m:j+m ) = J( i-m:i+m, j-m:j+m ) + Z * w;
        
        weight(i, j) =  sum( Z(:).*Z(:) );
        
    end
end

J = I - J / n^2;

end
