function [ sino, geom, secondSinoIndex] = mergeKvmvShortScans( sino1, sino2, geom1, geom2, totalAngle )
% [ sino, secondSinoIndex] = mergeOrthogonalShortScans( sino1, sino2, ratio )
% 
% Note: 
%       1. two sinogram must have exactly same geomety
%           can be done by function convertSinogramGeometry( )
%
%       2. secondSinoIndex - secondSinoIndex of the first projection that uses second sinogram 
% Meng Wu @ Stanford Univercity
% 2014.2


totalAngle = totalAngle * pi / 180;

fprintf( 'Merge two orthogonal short scan sinograms for %f degrees short scan. \n', totalAngle );

% rebin the second detector pixel to the same pixel width as the first one
sino2 = rebinFlatPanelDetector( sino2, geom2, geom1 );

% find the start and stop angles
middleAngle = ( max( geom1.betas) + min(geom2.betas) ) / 2;

startAngle = middleAngle - ( totalAngle - pi /2);
stopAngle = middleAngle + pi / 2;

geom = geom1;

% find the segment of sinograms
startIndex1 = find( geom1.betas > startAngle &  geom1.betas < middleAngle , 1 );
stopIndex1 = find( geom1.betas > startAngle &  geom1.betas < middleAngle , 1, 'last' ) ;

startIndex2 = find( geom2.betas > middleAngle &  geom2.betas <= stopAngle  , 1 );
stopIndex2 = find( geom2.betas > middleAngle &  geom2.betas <= stopAngle  , 1, 'last' );

geom.noViews1   = stopIndex1 - startIndex1 + 1;
geom.noViews2   = stopIndex2 - startIndex2 + 1;
geom.noViews    = geom.noViews1 + geom.noViews2;
geom.couchZ     = zeros( 1, geom.noViews );


secondSinoIndex = geom.noViews1 + 1;

geom.betas = zeros(1, geom.noViews);
geom.betas( 1:geom.noViews1 )     = geom1.betas( startIndex1:stopIndex1 );
geom.betas( geom.noViews1+1:end ) = geom2.betas( startIndex2:stopIndex2 );
geom.secondSinoIndex = secondSinoIndex;

% angular step size
geom.dbeta1 = abs( geom1.betas(stopIndex1) - geom1.betas(startIndex1) ) / ( stopIndex1 - startIndex1 );
geom.dbeta2 = abs( geom2.betas(stopIndex2) - geom2.betas(startIndex2) ) / ( stopIndex2 - startIndex2 );
geom.dbetaAverage = abs( geom.betas(end) - geom.betas(1) ) / ( geom.noViews - 1 );


sino = zeros( [geom.detSize(2), geom.detSize(1), geom.noViews ], 'single' );

if length( size( squeeze( sino ) ) ) == 3
    sino(:,:,1 : geom.noViews1)     = sino1(:,:,startIndex1:stopIndex1 );
    sino(:,:,geom.noViews1+1:end )  = sino2(:,:,startIndex2:stopIndex2);
else
    sino(:,1 : geom.noViews1)     = sino1(:,startIndex1:stopIndex1 );
    sino(:,geom.noViews1+1:end )  = sino2(:,startIndex2:stopIndex2);
end


fprintf('Done!\n\n\n')

end