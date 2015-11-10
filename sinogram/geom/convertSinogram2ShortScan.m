function [geomShort, sinoShort ]  = convertSinogram2ShortScan( geom, sino, offsetAngle )
% function [geomShort, sinoShort ]  = convertSinogram2ShortScan( geom, sino )
% note:
%       offsetAngle ( in degrees )
%
% Meng Wu @ Stanford Univercity
% 2014.2

if nargin < 3
    offsetAngle = 0;
end


fprintf( 'Convert full scna to short scan geometry with offest angle = %f degree. \n', offsetAngle );

offsetAngle = offsetAngle * pi / 180;
geomShort = geom;
% compute fan angles

miniangle = min( geom.betas );

geom.betas = geom.betas - miniangle;

if geom.flatPanel
    fanAngle = 2 * atan ( ( geom.detSpacing(1) * geom.detSize(1) * 0.55 ) / geom.SDD ) ;
else
    fanAngle = 2 * ( geom.detSpacing(1) * geom.detSize(1) * 0.55 ) / geom.SDD  ;
end

fprintf( '\tSystem fan angle = %2.3f degree. \n', 180 * fanAngle / pi );

% deal with speical case
if offsetAngle >=  pi - fanAngle
    
    fprintf('The offset angle is larger than 180 degrees.\n');
    if ndims( sino) == 3
        sinoTemp = sino;
        sino(:,:,1:end/2 )      = sinoTemp(:,:,end/2+1:end);
        sino(:,:,end/2+1:end)   = sinoTemp(:,:,1:end/2);
    else
        sinoTemp = sino;
        sino(:,1:end/2 )      = sinoTemp(:,end/2+1:end);
        sino(:,end/2+1:end)   = sinoTemp(:,1:end/2);
    end
    
    betasTemp = geom.betas;
    geom.betas(1:end/2 )    = betasTemp(end/2+1:end);
    geom.betas(end/2+1:end) = betasTemp(1:end/2) + pi * 2;
    
    couchZTemp = geom.couchZ;
    geom.couchZ(1:end/2 )    = couchZTemp(end/2+1:end);
    geom.couchZ(end/2+1:end) = couchZTemp(1:end/2);
    
end

% extract short scan data
viewStart = find( geom.betas >= offsetAngle, 1, 'first'  );
viewStop  = find( geom.betas <= offsetAngle + fanAngle + pi, 1, 'last'  ) + 1;

geomShort.betas     = geomShort.betas( viewStart:viewStop  ) ;
geomShort.couchZ    = geomShort.couchZ( viewStart:viewStop  );
geomShort.noViews   = viewStop - viewStart + 1;
geomShort.shortScan = 1;

if ndims( sino) == 3
    sinoShort = sino( :, :, viewStart:viewStop );
else
    sinoShort = sino( :, viewStart:viewStop );
end

fprintf('\tGeomerty information:\n');
fprintf('\tSDD: %.0fmm, SAD: %.0fmm\n', geomShort.SDD, geom.SAD);
fprintf('\tProjections: %i form %0.2f to %0.2f degrees\n', geomShort.noViews, min(geomShort.betas) * 180 / pi, max(geomShort.betas) * 180 / pi );
fprintf('\tDetector size: %ipx X %ipx\n',  geomShort.detSize(1), geomShort.detSize(end));
fprintf('\tPixel size: %5.3fmm X %2.3f mm  \n',  geomShort.detSpacing(1), geomShort.detSpacing(end) );
fprintf('\tReconstruction size: %i  X %i  X %i  \n',  geom.reconSize(1),  geom.reconSize(2), geom.reconSize(end) );
fprintf('\tReconstruction spacing: %2.3fmm X %2.3f mm X %2.3f mm  \n',  geom.reconSpacing(1), geom.reconSpacing(2), geom.reconSpacing(end) );

fprintf('Done!\n\n\n')

end