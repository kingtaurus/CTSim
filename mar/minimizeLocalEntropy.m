function [J,  weight] = minimizeLocalEntropy( I, C, n, bin, map )

m = (n - 1)/2;

weight = zeros( size(I) );
J = zeros( size(I) );

for i = n : size(I, 1) - n
    for j =  n : size( I, 2) - n

        X = I( i-m:i+m, j-m:j+m );
        Y = C( i-m:i+m, j-m:j+m );
        
        if  map(i,j) == 0 || max(Y(:)) < 2 * bin || max(Y(:)) > 100 * bin
            continue;
        end
        
        X = ( X - mean(X(:)) )/ bin / 256 + 0.5;
        Y = ( Y - mean(Y(:)) )/ bin / 256;
        Smin = inf;
        wmin = 0;
        
        for w = -0.05 : 0.01 : 0.05
            Z = X - Y * w;
            S = entropy( Z );
            if S < Smin
                Smin = S;
                wmin = w;
            end
        end
        
        J( i-m:i+m, j-m:j+m ) = J( i-m:i+m, j-m:j+m ) + Z * wmin;
        weight(i, j) = wmin;
        
    end
end

J = J / n^2 ;

end