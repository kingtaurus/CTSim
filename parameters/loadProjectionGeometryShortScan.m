function [ geom ] = loadProjectionGeometryShortScan( p, scanAngle, offsetAngle )
% Load short scan system geometry
%   inputs: 
%       p - configuration parameters
%       scanAnge - total short scan angle (in degree)
%       offsetAngle - offset of the starting angle (default 0)
%   output:
%       geom    - geometery parameters
%
% Copyright (c) 2010-2012 by Andreas Keil, Stanford University.
% Modified by Meng Wu at 2012.9

if nargin < 2
    scanAngle = 200;
end

if nargin < 3
    offsetAngle = 0;
end 

fprintf('Loading system geometry ...');

fprintf('Scan angle is %f degrees, offset angle is %f degree. \n', scanAngle, offsetAngle);

scanAngle = degtorad(scanAngle);
offsetAngle = degtorad(offsetAngle);


% origin coordinate of the reconstruction (world coordinate of the first voxel in mm)
geom.originRecon = -(p.Reconstruction.size-1)/2 .* p.Reconstruction.spacing;

% projection geometry
geom.SAD = p.Geometries.SAD;
geom.ADD = p.Geometries.ADD;
geom.SDD = geom.SAD + geom.ADD; % source-to-detector distance (mm)
geom.detSize     = p.Geometries.sizeDet;
geom.detSpacing  = p.Geometries.spacingDet;
geom.detOffset   = p.Geometries.offsetDet;
geom.flatPanel   = p.Geometries.flatPanel;
geom.noViews     = p.Geometries.noViews;
geom.couchZ      = zeros(1, geom.noViews, 'single');

geom.DQE         = p.Detector.detectorConversionEfficiency;
geom.reconSize   = p.Reconstruction.size;
geom.reconSpacing= p.Reconstruction.spacing;
geom.reconOffset = p.Reconstruction.offset;
geom.betas       = (0:geom.noViews-1)*scanAngle/(geom.noViews-1) + offsetAngle - scanAngle/2;

geom.shortScan      = 1;
geom.helicalScan    = 0;

printOutGeometryInfo( geom );
fprintf( 'done.\n\n');

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


if( pi + 2 * halfFanAngle > scanAngle )
    fprintf('Warning: scan angle is smaller then pi+fan. \n');
end

x = ( -(geom.reconSize(1)-1)/2:(geom.reconSize(1)-1)/2);
y = ( -(geom.reconSize(2)-1)/2:(geom.reconSize(2)-1)/2);
[xx, yy] = meshgrid( y, x);
geom.map = sqrt( xx.^2 + yy.^2 ) < max( max( geom.reconSize ), geom.FOV/geom.reconSpacing(1) ) /2;

end


function rad = degtorad(deg)

rad = pi * deg / 180;

end