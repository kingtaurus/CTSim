function dx = gradientRatioUpdata( gradientWLS, gradientPenalty, diff, forwardPathSeeking, dv, updateRate, totalNumUpdatePixels )

dx = zeros( size( gradientWLS ), 'single' );

% pixel with the same gradient direction
pixelSameDirection = gradientWLS .*  gradientPenalty > 0 ;

% normalize the gradiants
gradientWLS = gradientWLS / std(gradientWLS(:));
gradientPenalty = gradientPenalty / std(gradientPenalty(:));

if forwardPathSeeking
    lambda = gradientPenalty ./ ( abs( gradientWLS ) + 1e-3 );
    lambda( pixelSameDirection ) = lambda( pixelSameDirection ) * 100;
else
    lambda = gradientWLS ./ ( abs( gradientPenalty ) + 1e-3 );
    lambda( pixelSameDirection ) = lambda( pixelSameDirection ) * 100;
end

% make sure the update pixel are in the right direction
valid = abs( diff ) >  dv;
valid = valid & ( lambda .* diff < 0 ) ;

% rise up the percentage a little bit if only a few valid pixels
q = min( totalNumUpdatePixels * updateRate / sum(valid(:)), 0.5 );

% find the 1-q quantile on the scores
t = quantile( abs( lambda(valid(:)) ), 1 - q );

updatePixelMap       = ( abs( lambda ) > t ) & valid;
pixelUpdateFraction  = sum( updatePixelMap(:) ) / totalNumUpdatePixels;

% fixed step update in the gradient descent
dx( updatePixelMap ) = - dv * sign( diff( updatePixelMap ) );
dx = dx - dv * ( updateRate -  pixelUpdateFraction )  / mean( abs( diff(:) ) ) * diff ;

end