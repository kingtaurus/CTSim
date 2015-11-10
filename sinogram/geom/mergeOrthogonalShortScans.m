function [ sino, secondSinoIndex] = mergeOrthogonalShortScans( sino1, sino2, ratio, geom )
% [ sino, secondSinoIndex] = mergeOrthogonalShortScans( sino1, sino2, ratio )
% 
% Note: 
%       1. two sinogram must have exactly same geomety
%           can be done by function convertSinogramGeometry( )
%
%       2. secondSinoIndex - secondSinoIndex of the first projection that uses second sinogram 
% Meng Wu @ Stanford Univercity
% 2014.2



fprintf( 'Merge two orthogonal short scan sinograms with ratio = %f. \n', ratio );

sino = sino1;
if length( size( squeeze( sino ) ) ) == 3
    noView = size( sino, 3);
    secondSinoIndex = round( ratio * noView ) + 1;
    sino(:,:,secondSinoIndex:end ) = sino2(:,:,secondSinoIndex : end);
else
    noView = size( sino, 2);
    secondSinoIndex = round( ratio * noView ) + 1;
    sino(:,secondSinoIndex:end ) = sino2(:,secondSinoIndex : end );
end


if nargin > 3
    minScanAngle = max( ratio, 1 - ratio ) * abs( geom.betas(end) - geom.betas(1) - pi / 2 );
    minScanAngle = minScanAngle * 180 / pi;
    minDuration = minScanAngle / 6 ;
    
    fprintf( '\tThe minimum scan anlge is %3.1f degrees, minimum time is %2.1f sec. \n', minScanAngle, minDuration );
    
end

geom.noViews1   = secondSinoIndex - 1;
geom.noViews2   = noView - secondSinoIndex + 1;

% angular step size
geom.dbetaAverage = abs( geom.betas(end) - geom.betas(1) ) / ( geom.noViews - 1 );
geom.dbeta1 = geom.dbetaAverage;
geom.dbeta2 = geom.dbetaAverage;


fprintf('Done!\n\n\n')

end