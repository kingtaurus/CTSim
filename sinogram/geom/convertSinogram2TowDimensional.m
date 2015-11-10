function [geom2d, sino2d ]  = convertSinogram2TowDimensional( geom, sino )
% function [geom2d, sino2d ]  = convertSinogram2TowDimensiontail( geom, sino )
% note:
%       offsetAngle ( in degrees )
%
% Meng Wu @ Stanford Univercity
% 2014.2


fprintf( 'Convert 3d scan data to 2d scan. \n' );


sino2d = squeeze( sino(round(end/2),:,:) );
geom2d = geom;
geom2d.reconSize 	= geom2d.reconSize(1:2);
geom2d.reconSpacing = geom2d.reconSpacing(1:2);
geom2d.reconOffset = geom2d.reconOffset(1:2);
geom2d.detSize     = geom2d.detSize(1);
geom2d.detSpacing  = geom2d.detSpacing(1);
geom2d.detOffset   = geom2d.detOffset(1);


fprintf('\tGeomerty information:\n');
fprintf('\tSDD: %.0fmm, SAD: %.0fmm\n', geom2d.SDD, geom.SAD);
fprintf('\tProjections: %i form %0.2f to %0.2f degrees\n', geom2d.noViews, min(geom2d.betas) * 180 / pi, max(geom2d.betas) * 180 / pi );
fprintf('\tDetector size: %ipx X %ipx\n',  geom2d.detSize(1), geom2d.detSize(end));
fprintf('\tPixel size: %5.3fmm X %2.3f mm  \n',  geom2d.detSpacing(1), geom2d.detSpacing(end) );
fprintf('\tReconstruction size: %i  X %i  X %i  \n',  geom.reconSize(1),  geom.reconSize(2), geom.reconSize(end) );
fprintf('\tReconstruction spacing: %2.3fmm X %2.3f mm X %2.3f mm  \n',  geom.reconSpacing(1), geom.reconSpacing(2), geom.reconSpacing(end) );

fprintf('Done!\n\n\n')

end