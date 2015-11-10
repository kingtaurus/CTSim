function beta = estimateBetaKKT(   gradientWLS, gradientPenalty, validPixelsAll, beta0, forwardPathSeeking )
% using the KKT condition to estimate the tuning parameter value

gradientRatio   =  - gradientWLS ./ gradientPenalty;
valid           = validPixelsAll & ~isnan( gradientRatio ) ...
    & ~isinf( gradientRatio ) & abs( gradientPenalty ) > 1e-8;

beta =  median( gradientRatio(valid(:)) );

% get a little extreme
if forwardPathSeeking
    beta = max( 1.05 * beta, beta0 );
else
    beta = min( 0.95 * beta, beta0 );
end

end
