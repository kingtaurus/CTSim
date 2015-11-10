function [ geomKeV, geomMeV] = loadProjectionGeometryKvmv( p )
% Load system geometry
%
% Copyright (c) 2010-2012 by Andreas Keil, Stanford University.
% Modified by Meng Wu at 2012.9

fprintf('Loading system geometry ...');

% origin coordinate of the reconstruction (world coordinate of the first voxel in mm)
geomKeV.originRecon = -(p.Reconstruction.size-1)/2 .* p.Reconstruction.spacing;
geomMeV.originRecon = geomKeV.originRecon;

% projection geometry
%geomKeV.allBetas    = (0:p.Geometries.keV.noViews-1) * 2*pi/p.Geometries.keV.noViews;
geomKeV.SAD         = p.Geometries.keV.SAD;
geomKeV.ADD         = p.Geometries.keV.ADD;
geomKeV.SDD         = geomKeV.SAD + geomKeV.ADD; % source-to-detector distance (mm)
geomKeV.detSize     = p.Geometries.keV.sizeDet;
geomKeV.detSpacing  = p.Geometries.keV.spacingDet;
geomKeV.detOffset   = p.Geometries.keV.offsetDet;
geomKeV.noViews     = p.Geometries.keV.noViews;
%geomKeV.DQE         = p.Detector.detectorConversionEfficiencyKeV;
geomKeV.reconSize   = p.Reconstruction.size;
geomKeV.reconSpacing= p.Reconstruction.spacing;
geomKeV.reconOffset = p.Reconstruction.offset;
geomKeV.betas       = (0 : geomKeV.noViews-1) * 2 * pi / geomKeV.noViews;
geomKeV.flatPanel   = 1;
geomKeV.couchZ      = zeros(1, geomKeV.noViews);
% metal segmentation methods
geomKeV.sinoBasedSegm =  p.Geometries.keV.sinoBasedSegm;

geomKeV.shortScan      = 0;
geomKeV.helicalScan    = 0;


%geomMeV.allBetas    = (0:p.Geometries.MeV.noViews-1) * 2*pi/p.Geometries.MeV.noViews;
geomMeV.SAD         = p.Geometries.MeV.SAD;
geomMeV.ADD         = p.Geometries.MeV.ADD;
geomMeV.SDD         = geomMeV.SAD + geomMeV.ADD; % source-to-detector distance (mm)
geomMeV.detSize     = p.Geometries.MeV.sizeDet;
geomMeV.detSpacing  = p.Geometries.MeV.spacingDet;
geomMeV.detOffset   = p.Geometries.MeV.offsetDet;
geomMeV.noViews     = p.Geometries.MeV.noViews;
%geomMeV.DQE         = p.Detector.detectorConversionEfficiencyMeV;
geomMeV.reconSize   = p.Reconstruction.size;
geomMeV.reconSpacing= p.Reconstruction.spacing;
geomMeV.reconOffset = p.Reconstruction.offset;
geomMeV.betas       = (0 : geomMeV.noViews-1) * 2 * pi / geomMeV.noViews;
geomMeV.flatPanel   = 1;
geomMeV.couchZ      = zeros(1, geomMeV.noViews);

geomMeV.shortScan      = 0;
geomMeV.helicalScan    = 0;


fprintf('\tKV Geomerty information:\n');
fprintf('\tSDD: %.0fmm, SAD: %.0fmm\n', geomKeV.SDD, geomKeV.SAD);
fprintf('\tProjections: %i form %0.2f to %0.2f degrees\n', geomKeV.noViews, min(geomKeV.betas) * 180 / pi, max(geomKeV.betas) * 180 / pi );
fprintf('\tDetector size: %ipx X %ipx\n',  geomKeV.detSize(1), geomKeV.detSize(end));
fprintf('\tPixel size: %5.3fmm X %2.3f mm  \n',  geomKeV.detSpacing(1), geomKeV.detSpacing(end) );
fprintf('\tReconstruction size: %i  X %i  X %i  \n',  geomKeV.reconSize(1),  geomKeV.reconSize(2), geomKeV.reconSize(end) );
fprintf('\tReconstruction spacing: %2.3fmm X %2.3f mm X %2.3f mm  \n',  geomKeV.reconSpacing(1), geomKeV.reconSpacing(2), geomKeV.reconSpacing(end) );


fprintf('\n\tMV Geomerty information:\n');
fprintf('\tSDD: %.0fmm, SAD: %.0fmm\n', geomMeV.SDD, geomMeV.SAD);
fprintf('\tProjections: %i form %0.2f to %0.2f degrees\n', geomMeV.noViews, min(geomMeV.betas) * 180 / pi, max(geomMeV.betas) * 180 / pi );
fprintf('\tDetector size: %ipx X %ipx\n',  geomMeV.detSize(1), geomMeV.detSize(end));
fprintf('\tPixel size: %5.3fmm X %2.3f mm  \n',  geomMeV.detSpacing(1), geomMeV.detSpacing(end) );
fprintf('\tReconstruction size: %i  X %i  X %i  \n',  geomMeV.reconSize(1),  geomMeV.reconSize(2), geomMeV.reconSize(end) );
fprintf('\tReconstruction spacing: %2.3fmm X %2.3f mm X %2.3f mm  \n',  geomMeV.reconSpacing(1), geomMeV.reconSpacing(2), geomMeV.reconSpacing(end) );



fprintf( 'done.\n\n');



end