load 'temp.mat'

% Load simulation parameters and datas
[phan, map] = loadXCATPhantom(p);

phan.offset = [0 90 50]

turns = 5;
pitch = 0.5;

[ geom ] = loadProjectionGeometryHelicalCT( p, turns, pitch );

sinosDirKeV = [sinosDirKeV sprintf( '%i-turns-%i-pitch-%i-shift', turns,...
    round( pitch * 100 ) , round( phan.offset(2) )), '/'];

spectrum = loadSpectraCT(p, geom, 2e5);

[sinoRaw, tubeCurrentProfile] = simulateCTRawData( phan, geom, spectrum, sinosDirKeV );

sinoAtt = processCTRawData( sinoRaw, spectrum, tubeCurrentProfile);


%%

mu = 0.2;
edgeWidth = 6;
extendWidth = 128;
reconFunc = @(sino, geom )reconHelical( sino, geom, 'hamming', 'wfdk', 'short' );

%%
img_no_truncation_crr = reconFunc( sinoAtt, geom );


sinoquick( sinoAtt );
imquick( img_no_truncation_crr, map.windowAtt);

%%

[ sinoAtt_ellp, geom_ellp ] = truncationCorrectionEllipicalPatch( sinoAtt, geom, mu, edgeWidth, 1, extendWidth );
img_truncation_crr_ellp = reconFunc( sinoAtt_ellp, geom_ellp );

sinoquick( sinoAtt_ellp );
imquick( img_truncation_crr_ellp, map.windowAtt);


%%

[ sinoAtt_wcf, geom_wcf ] = truncationCorrectionWaterCylinderFitting( sinoAtt, geom, mu, edgeWidth, extendWidth );
img_truncation_crr_wcf = reconFunc( sinoAtt_wcf, geom_wcf);

sinoquick( sinoAtt_wcf );
imquick( img_truncation_crr_wcf, map.windowAtt);

%%
[ sinoAtt_copy, geom_copy ] = truncationCorrectionCopyPatch( sinoAtt, geom, mu, edgeWidth, extendWidth );
img_truncation_crr_copy = reconFunc( sinoAtt_copy, geom_copy);

sinoquick( sinoAtt_copy );
imquick( img_truncation_crr_copy, map.windowAtt);


%%

[ sinoAtt_wcftc, geom_wcftc ] = truncationCorrectionWaterCylinderFittingTotalConsistency( sinoAtt, geom, mu, edgeWidth, extendWidth );
img_truncation_crr_wcftc = reconFunc( sinoAtt_wcftc, geom_wcftc);

sinoquick( sinoAtt_wcftc );
imquick( img_truncation_crr_wcftc, map.windowAtt);

%% 
[ sinoAtt_wcfs, geom_wcfs ] = truncationCorrectionWaterCylinderFittingSmoothConstrain( sinoAtt, geom, mu, edgeWidth, extendWidth );
img_truncation_crr_wcf_smooth = reconFunc( sinoAtt_wcfs, geom_wcfs);

sinoquick( sinoAtt_wcfs );
imquick( img_truncation_crr_wcf_smooth, map.windowAtt);


%%

[ sinoAtt_2pass, geom_2pass ] = truncationCorrectionTwoPass(  sinoAtt, geom, mu, edgeWidth, extendWidth, reconFunc );
img_truncation_crr_2pass = reconFunc( sinoAtt_2pass, geom_2pass);

sinoquick( sinoAtt_2pass );
imquick( img_truncation_crr_2pass, map.windowAtt);


%% Simulate reference detector
geomRef = geom;
geomRef.detSize = [ geom.detSize(1) + 2 * extendWidth,  2 ]; 

spectrumRef = loadSpectraCT(p, geomRef, 2e5);

[sinoRawRef, tubeCurrentProfile] = simulateCTRawData( phan, geomRef, spectrumRef, [ sinosDirKeV 'ref/'] );

sinoRef = processCTRawData( sinoRawRef, spectrumRef, tubeCurrentProfile);

%% combine regular sinogram with extended reference detector 

sinoAtt2 = zeros( [geom.detSize(2), geom.detSize(1) + 2 * extendWidth,  geom.noViews], 'single' );

sinoAtt2( :, extendWidth+1:end-extendWidth , : ) = sinoAtt;
sinoAtt2( end/2:end/2+1, :, : ) = sinoRef;

figure; imdisp( sinoAtt2(:,:,3000 )' );

geom2 = geomRef;
geom2.detSize(2) = geom.detSize(2);

%% reconstruction using reference detector 
img_ref = reconFunc( sinoRef, geomRef );
sinoquick( sinoRef );
imquick( img_ref, map.windowAtt);

%%
sinoRefLeft = squeeze( mean( sinoRef(:, 1 : extendWidth, : ), 1) );
sinoRefRight = squeeze( mean( sinoRef(:, end - extendWidth + 1 : end, : ), 1) );

[ sinoAtt_ref, geom_ref ] = truncationCorrectionReferenceDetector(  sinoAtt, geom, sinoRefLeft, sinoRefRight, edgeWidth, extendWidth );
img_truncation_crr_ref = reconFunc( sinoAtt_ref, geom_ref);

sinoquick( sinoAtt_ref );
imquick( img_truncation_crr_ref, map.windowAtt);

%%

[ geom ] = loadProjectionGeometryHelicalCT( p, turns, pitch );

sinosDirKeV = [sinosDirKeV sprintf( '%i-turns-%i-pitch-%i-shift', turns,...
    round( pitch * 100 ) , round( phan.offset(2) )), '/'];

spectrum = loadSpectraCT(p, geom, 2e5);

[sinoRaw, tubeCurrentProfile] = simulateCTRawData( phan, geom, spectrum, sinosDirKeV );

sinoAtt = processCTRawData( sinoRaw, spectrum, tubeCurrentProfile);

%% Simulate reference detector
[ geomFlat ] = loadProjectionGeometryCT( p );

geomFlat.flatPanel = 1;

spectrumFlat = loadSpectraCT(p, geomFlat, 2e5);

[sinoRawFlat, tubeCurrentProfile] = simulateCTRawData( phan, geomFlat, spectrumFlat, [ sinosDirKeV 'flat/'] );

sinoFlat = processCTRawData( sinoRawFlat, spectrumFlat, tubeCurrentProfile);


%% ATRACT algorithm
geomFlat.reconSize(3) = 8;
geomFlat.shortScan = 0;
img_atract = reconATRACT1( sinoFlat, geomFlat );
sinoquick( sinoFlat );
imquick( img_atract, [0 80]);
