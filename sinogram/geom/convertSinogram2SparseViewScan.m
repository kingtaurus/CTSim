function [geomSparse, sinoSparse, weightsSparse ]  = convertSinogram2SparseViewScan( geom, sino, rate, weights )
% function [geomSparse, sinoSparse ]  = convertSinogram2SparseViewScan( geom, sino, rate )
% note:
%       rate - downsample rate
%
% Meng Wu @ Stanford Univercity
% 2014.2


fprintf( 'Convert to sparse view scan geometry with downsample rate = %f. \n', rate );

geomSparse = geom ;
geomSparse.betas     = geomSparse.betas( 1:rate:end );
geomSparse.couchZ    = geomSparse.couchZ( 1:rate:end );
geomSparse.noViews   = length( geomSparse.betas );

if isfield( geom, 'SADPerProjection' )
    geomSparse.SADPerProjection = geom.SADPerProjection( 1:rate:end );
end

if isfield( geom, 'SDDPerProjection' )
    geomSparse.SDDPerProjection = geom.SDDPerProjection( 1:rate:end );
end

if isfield( geom, 'ADDPerProjection' )
    geomSparse.ADDPerProjection = geom.ADDPerProjection( 1:rate:end );
end

if isfield( geom, 'detOffsetsPerProjection' )
    geomSparse.detOffsetsPerProjection = geom.detOffsetsPerProjection( :, 1:rate:end );
end

if isfield( geom, 'cameraPositions' )
    geomSparse.cameraPositions = geom.cameraPositions( :, 1:rate:end );
end

if isfield( geom, 'PMat' )
    geomSparse.PMat = geom.PMat( :, :, 1:rate:end );
end


if ndims( sino) == 3
    sinoSparse = sino( :, :, 1:rate:end );
else
    sinoSparse = sino( :, 1:rate:end );
end

if nargin == 4
    if ndims( sino) == 3
        weightsSparse = weights( :, :, 1:rate:end );
    else
        weightsSparse = weights( :, 1:rate:end );
    end
else
    weightsSparse = 0;
end


fprintf('\tGeomerty information:\n');
fprintf('\tSDD: %.0fmm, SAD: %.0fmm\n', geomSparse.SDD, geom.SAD);
fprintf('\tProjections: %i form %0.2f to %0.2f degrees\n', geomSparse.noViews, min(geomSparse.betas) * 180 / pi, max(geomSparse.betas) * 180 / pi );
fprintf('\tDetector size: %ipx X %ipx\n',  geomSparse.detSize(1), geomSparse.detSize(end));
fprintf('\tPixel size: %5.3fmm X %2.3f mm  \n',  geomSparse.detSpacing(1), geomSparse.detSpacing(end) );
fprintf('\tReconstruction size: %i  X %i  X %i  \n',  geom.reconSize(1),  geom.reconSize(2), geom.reconSize(end) );
fprintf('\tReconstruction spacing: %2.3fmm X %2.3f mm X %2.3f mm  \n',  geom.reconSpacing(1), geom.reconSpacing(2), geom.reconSpacing(end) );

fprintf('Done!\n\n\n')

end