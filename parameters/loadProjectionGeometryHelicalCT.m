function [ geom ] = loadProjectionGeometryHelicalCT( p, noTurns, pitch )
% Load helical CT system geometry
%   inputs:
%       p - configuration parameters
%       noTurns - number of hilical turns
%       pitch   - helical pitch ( couch movement per 360 scan / detector height)
%   output:
%       geom    - geometery parameters
%
% Meng Wu at Stanford University
% 2014.3

fprintf('Loading system geometry ...\n');

% origin coordinate of the reconstruction (world coordinate of the first voxel in mm)
geom.originRecon = -(p.Reconstruction.size-1)/2 .* p.Reconstruction.spacing;

% projection geometry
geom.SAD        = p.Geometries.SAD;
geom.ADD        = p.Geometries.ADD;
geom.SDD        = geom.SAD + geom.ADD; % source-to-detector distance (mm)
geom.detSize     = p.Geometries.sizeDet;
geom.detSpacing  = p.Geometries.spacingDet;
geom.detOffset   = p.Geometries.offsetDet;

geom.noTurns        = noTurns;
geom.noViewsTurn    = p.Geometries.noViews;
geom.noViews        = p.Geometries.noViews * noTurns + 1;

geom.pitch          = pitch;
geom.couchSpacing   = pitch * geom.detSize(2) * geom.detSpacing(2) / geom.noViewsTurn;
geom.couchZ         = ( ( 1: geom.noViews ) -  ( geom.noViews + 1 ) / 2 ) * geom.couchSpacing;

geom.flatPanel      = p.Geometries.flatPanel;
geom.DQE            = p.Detector.detectorConversionEfficiency;

geom.reconSize      = p.Reconstruction.size;
geom.reconSpacing   = p.Reconstruction.spacing;
geom.reconOffset    = p.Reconstruction.offset;
geom.betas          = ( ( 1: geom.noViews ) -  ( geom.noViews + 1 ) / 2 ) * 2 * pi / geom.noViewsTurn + pi;

geom.shortScan      = 1;
geom.helicalScan    = 1;

if geom.flatPanel
    halfFanAngle   = atan( ( geom.detSize(1) * geom.detSpacing(1) ) / 2 / geom.SDD );
    geom.FOV       = 2 * sin( halfFanAngle ) * geom.SAD;
else
    halfFanAngle   = ( geom.detSize(1) * geom.detSpacing(1) ) / 2 / geom.SDD ;
    geom.FOV       = 2 * sin( halfFanAngle ) * geom.SAD;
end

x = ( -(geom.reconSize(1)-1)/2:(geom.reconSize(1)-1)/2);
y = ( -(geom.reconSize(2)-1)/2:(geom.reconSize(2)-1)/2);

[xx, yy] = meshgrid( y, x);
geom.map = sqrt( xx.^2 + yy.^2 ) < max( max( geom.reconSize ), geom.FOV/geom.reconSpacing(1) ) /2;

geom.detPSF     = p.Detector.pointSpreadFunctionFWHM ;
geom.detNPS     = p.Detector.noisePowerSpectrum;
geom.focalPSF   = p.Spectra.focalSpotSize ;


printOutGeometryInfo( geom );

fprintf( 'done.\n\n');

end