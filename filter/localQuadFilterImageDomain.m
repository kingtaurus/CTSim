function imgOut = localQuadFilterImageDomain( imgIn, imgVar, w, beta )

fprintf('Local quadratic programming filter in image domain with beta = %2.2f ... \n', beta );


imgOut = zeros( size(imgIn), class(imgIn) );

for iz = 1 : size( imgIn, 3 )
    
    imgOut(:,:,iz) = LQPFilter( squeeze( imgIn(:,:,iz) ), sqrt( squeeze( imgVar(:,:,iz) ) ), w, 1 );
    
end

fprintf('\t done.\n');

end

function [uhat] = LQPFilter(v, s, wn, beta)
% get the size of the image
options = optimset('quadprog');
options = optimset(options,'Algorithm', 'interior-point-convex',...
    'Display','off', 'MaxIter', 10);

v( isnan(v) ) = -1000;

n = length(v);

u = v;
%u = adaptiveBilateralFilter( v, 2, 1, s );

% iterate over the pixels
ut = zeros(n);
%h = waitbar(0,'message') ;
for x = 1:n
    for y = 1:n
        % compute the search window
        iMin = max(x-wn,1);
        iMax = min(x+wn,n);
        jMin = max(y-wn,1);
        jMax = min(y+wn,n);
        
        V = v(iMin:iMax , jMin:jMax);
        U = u(iMin:iMax , jMin:jMax) - u(x,y);
        S = s(iMin:iMax , jMin:jMax);

        if s( x, y) == 0 || sum( isnan( U(:) ) ) || sum( isnan( V(:) ) )
            continue;
        end
        
        m = length(V(:));
        z = U(:);
        
        %w = exp( - z.^2 / ( 2 * s(x,y)^2 ) );
        %z = z - 0.8 * sum( z.* w ) / sum( w );
        %z = z .* max(abs(z) - s/2, 0) ./ ( abs(z) + 1);
        
        H = max( z*z' - s(x,y)^2, 0 ) + diag( S(:).^2 );
        w = quadprog( double(H ),[],[],[],ones(1,m),1,[],[],[],options );
        
        w = w / sum(w);
        
        % compute the weighted average over the search neighborhood
        ut(x,y) = sum(w(:) .* V(:));
        
    end
   % waitbar(x / n, h);
end
%close(h);
uhat = ut;

end