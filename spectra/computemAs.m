    load 'temp.mat';

geom = loadProjectionGeometryCT( p );

n0 = 1e6;
spectrum = loadSpectraCT(p, geom, n0);

tubeCurrentProfileFilename = 'tubeCurrentProfile_v4-XCAT-lung_kVp119.0_FLUX1e+06_BOWT1_ACE1.00_EID1_NOISE1_sim2.raw';

fileID = fopen(tubeCurrentProfileFilename, 'r');
tubeCurrentProfile = fread(fileID, Inf, 'double');
fclose(fileID);

Q = 2.6750e+06;

%%
fprintf('Simulation 6 using %s \n', tubeCurrentProfileFilename );
fprintf('Photons per pixels \t%.3g\n', n0 );
fprintf( 'Pixel size: \t\t\t%.1f x %.1f \n', geom.detSpacing(1), geom.detSpacing(2) );
fprintf( 'Simulation q / mm2 @1m: \t%.3g \n' , spectrum.photonsPerMm2At1m ); 
fprintf( '1 mAs  q / mm2 @1m: \t\t%.3g \n' , Q ); 
fprintf( 'Number of projections: \t\t%i \n', length( tubeCurrentProfile ) ); 
nrot = length( tubeCurrentProfile ) / geom.noViews;
fprintf( 'Number of rotation: \t\t%.1f \n', nrot );
fprintf( 'Average mAs per projection: \t%f \n',  mean( tubeCurrentProfile ) * spectrum.photonsPerMm2At1m / Q );
fprintf( 'Average mAs per rotation: \t%f \n',  sum( tubeCurrentProfile ) * spectrum.photonsPerMm2At1m / Q / nrot );
fprintf( 'Total mAs: \t\t\t%f \n\n\n',  sum( tubeCurrentProfile ) * spectrum.photonsPerMm2At1m / Q );
