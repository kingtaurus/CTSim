function [ geom ] = loadProjectionGeometryCT( p )
% Load system geometry
%   inputs: 
%       p - configuration parameters
%   output:
%       geom    - geometery parameters
%
% Copyright (c) 2010-2012 by Andreas Keil, Stanford University.
% Modified by Meng Wu at 2012.9

fprintf('Loading system geometry ...\n');

% origin coordinate of the reconstruction (world coordinate of the first voxel in mm)
geom.originRecon = -(p.Reconstruction.size-1)/2 .* p.Reconstruction.spacing;

% projection geometry
geom.SAD = p.Geometries.SAD;
geom.ADD = p.Geometries.ADD;
geom.SDD = geom.SAD + geom.ADD; % source-to-detector distance (mm)
geom.detSize     = p.Geometries.sizeDet;
geom.detSpacing  = p.Geometries.spacingDet;
geom.detOffset   = p.Geometries.offsetDet;
geom.noViews     = p.Geometries.noViews;
geom.flatPanel   = p.Geometries.flatPanel;
geom.DQE         = p.Detector.detectorConversionEfficiency;
geom.couchZ      = zeros(1, geom.noViews);

geom.reconSize   = p.Reconstruction.size;
geom.reconSpacing= p.Reconstruction.spacing;
geom.reconOffset = p.Reconstruction.offset;
%geom.betas       = (1 : geom.noViews) * 2 * pi / geom.noViews + pi ;
%geom.betas       = ( ( 1: geom.noViews ) -  ( geom.noViews + 1 ) / 2 ) * 2 * pi / geom.noViews ;
 geom.betas       = ( ( 1: geom.noViews ) -   geom.noViews / 2 ) * 2 * pi / geom.noViews ;
% geom.betas       = ( ( 1: geom.noViews ) -   geom.noViews / 2 - 1./25) * 2 * pi / geom.noViews;
geom.shortScan      = 0;
geom.helicalScan    = 0;

geom.detPSF = p.Detector.pointSpreadFunctionFWHM ;
geom.detNPS = p.Detector.noisePowerSpectrum;
geom.focalPSF = p.Spectra.focalSpotSize ;

if geom.flatPanel
    halfFanAngle = atan( ( geom.detSize(1) * geom.detSpacing(1) ) / 2 / geom.SDD );
    geom.FOV = 2 * sin( halfFanAngle ) * geom.SAD;
else
    halfFanAngle = ( geom.detSize(1) * geom.detSpacing(1) ) / 2 / geom.SDD ;
    geom.FOV = 2 * sin( halfFanAngle ) * geom.SAD;
end

x = ( -(geom.reconSize(1)-1)/2:(geom.reconSize(1)-1)/2);
y = ( -(geom.reconSize(2)-1)/2:(geom.reconSize(2)-1)/2);
[xx, yy] = meshgrid( y, x);
geom.map = sqrt( xx.^2 + yy.^2 ) < max( max( geom.reconSize ), geom.FOV/geom.reconSpacing(1) ) /2;


printOutGeometryInfo( geom );
fprintf( 'done.\n\n');



end