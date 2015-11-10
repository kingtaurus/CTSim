function helicalWeights = computeHelicalWeights3D( geom, noViewsPerSlice )
% function redundantWeights = computeHelicalWeights2D( gammas, betas, iview, nv, segmentLength )
%
%
% Meng Wu, 2014.7


% gat a sample sequences of rotation angle
betas = geom.betas( 1:noViewsPerSlice );
% center the roation
betas = betas - min( betas  );
% total the rotation angle
Delta = abs( betas(end) -  betas(1) );

angle2height = geom.pitch * geom.detSpacing(2) * geom.detSize(2) / ( 2 * pi );

u = ( ( -(geom.detSize(1)-1)/2:(geom.detSize(1)-1)/2) + geom.detOffset(1) ) * geom.detSpacing(1);
v = ( ( -(geom.detSize(2)-1)/2:(geom.detSize(2)-1)/2) + geom.detOffset(2) ) * geom.detSpacing(2);
[uu, vv ] = meshgrid( u, v);

% compute the fan angle of each detector column
if geom.flatPanel
    gammas   = atan( - u / geom.SDD );
    tanThetas   = vv ./ sqrt( uu.^2 + geom.SDD^2 ) ;
else
    gammas   = - u / geom.SDD ;
    tanThetas   = vv ./ geom.SDD ;
end

tanThetaMax = 1.1 * max( abs( tanThetas(:) ) );
k = 0.5  * pi / ( 2 * tanThetaMax );


helicalWeights = zeros( [ geom.detSize(2), geom.detSize(1), noViewsPerSlice], 'single' );

% for each projection angle
for iview = 1 : noViewsPerSlice
    
    if noViewsPerSlice < geom.noViewsTurn
        fanWeights = shortScanSliverWeight( gammas, betas(iview),  Delta - pi );
    else
        fanWeights = 0.5 * ones( size( gammas ) );
    end
    
    lambda =  betas( iview ) - Delta / 2 ;
    viewWeights = zeros( [ geom.detSize(2), geom.detSize(1)], 'single' );
    
    % for each dector column
    for i = 1:length(gammas)
        
        
        lambdaConj1 = lambda + pi - 2 * gammas(i);
        lambdaConj2 = lambda - pi - 2 * gammas(i);
        
        tanThetaColumn = tanThetas(:, i);
        
        
        if abs(lambdaConj1  ) < abs( lambdaConj2  )
            lambdaConj = lambdaConj1;
        else
            lambdaConj = lambdaConj2;
        end
        
        if abs( lambdaConj ) > Delta / 2 && noViewsPerSlice < geom.noViewsTurn
            viewWeights(:, i) = fanWeights(i);
            continue;
        end
        
        src2pixelDistance = geom.SDD - lambda * angle2height ./ tanThetaColumn ;
        
        src2pixelDistance( src2pixelDistance < 10 ) = 10;
        tanThetaConj = ( lambdaConj * angle2height ) ./  src2pixelDistance ;
        
        tanThetaConj( abs( tanThetaConj ) > tanThetaMax ) = tanThetaMax;
        rayWeights = cos( k * abs( tanThetaColumn ) );
        conjWeights = cos( k * abs( tanThetaConj  ) );
        
        rayWeights = rayWeights * fanWeights(i);
        conjWeights = conjWeights * ( 1 - fanWeights(i) );
        viewWeights(:, i) =  rayWeights  ./ ( rayWeights + conjWeights) ;
        
    end
    
    G = fspecial('gaussian',[1 8], 2);
    viewWeights = imfilter(viewWeights,G,'same');
    
    helicalWeights(:,:,iview) = viewWeights ;
    
    %imdisp( viewWeights' , [0 1])
end

end