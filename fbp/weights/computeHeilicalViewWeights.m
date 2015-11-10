function viewWeights = computeHeilicalViewWeights( gammas, detTanThetas, betas, iview, SDD, angle2height, segmentLength  )
%
%
% Meng Wu, 2014.3

Delta    = abs(betas(end) - betas(1));
beta     = betas(iview) - min(betas);
tanThetaMax =  1.1 * max( abs( detTanThetas ) );


if strcmpi(segmentLength, 'short')
    redundantWeights = shortScanSliverWeight( gammas, beta,  Delta - pi );
else
    redundantWeights = 0.5 * ones( size( gammas ) );
end


k = 0.5  * pi / ( 2 * tanThetaMax );

viewWeights = zeros( [length(detTanThetas), length(gammas)], 'single' );

lambda = beta - Delta / 2 ;

for i = 1:length(gammas)
    
    lambdaConj1 = lambda + pi - 2 * gammas(i);
    lambdaConj2 = lambda - pi - 2 * gammas(i);
    
    if abs(lambdaConj1  ) < abs( lambdaConj2  )
        lambdaConj = lambdaConj1;
    else
        lambdaConj = lambdaConj2;
    end
    
    if abs( lambdaConj ) > Delta / 2 
        viewWeights(:, i) = redundantWeights(i);
        continue;
    end
    
    src2pixelDistance = SDD - lambda * angle2height ./ detTanThetas ;
    tanThetaConj = ( lambdaConj * angle2height ) ./  src2pixelDistance ;
    
    tanThetaConj( abs( tanThetaConj ) > tanThetaMax ) = tanThetaMax;
    rayWeights = cos( k * abs( detTanThetas ) );
    conjWeights = cos( k * abs( tanThetaConj  ) );
    
    rayWeights = rayWeights * redundantWeights(i);
    conjWeights = conjWeights * ( 1 - redundantWeights(i) );
    viewWeights(:, i) =  rayWeights  ./ ( rayWeights + conjWeights + 1E-6 ) ;

end


viewWeights = medfilt2(viewWeights,[5 5]);

G = fspecial('gaussian',[7 7], 2);
viewWeights = imfilter(viewWeights,G,'replicate');


end