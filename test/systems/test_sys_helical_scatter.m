%% test_sys_helical_scatter.m
% This test script is mainly for testing the helical system geometry, 
% helical reconstruction, and scatter simulation.
%
% Meng Wu at Stanford University
% 2013 - 2014

%% Load simulation parameters and datas

load 'temp.mat';

% load phantom
[phan, map] = loadXCATPhantom(p);

turns = 6;
pitch = 31/32;
sinosDirKeV = [sinosDirKeV sprintf( '%i-turns-%i-pitch', turns, round( pitch * 100 ) ), '/'];

% Geometry info
geom = loadProjectionGeometryHelicalCT( p, turns, pitch );

% Specturm info
spectrum = loadSpectraCT(p, geom, 1e6);

% Compute ground truth
[imgGtAtt, imgGtHu ] = computeGroundTruth(phan, spectrum);


%% Simulate the raw data

% sinoRaw is the data set without any scatter
[sinoRaw, tubeCurrentProfile] = simulateCTRawData( phan, geom, spectrum, sinosDirKeV );

%% Load MC scatter simulation results

load( 'E:\MATLAB\CTData\ScatterSimulation\helical\scatter_helical_bowtie_pitch1_32rows.mat' );

% sinoRaw_scatter is the data set after adding the simulated scatter, no
% algorithm correction at all (anti-scatter grid is used)

[sinoRaw_scatter, sinoPrime, sinoScatter] = combineScatterSimulation( sinoRaw, primary_bowtie, scatter_bowtie, 10  );
% sinoRaw_scatter_corr is the data set after adding the simulated scatter,
% no algorithm correction at all (anti-scatter grid is used) 

sinoRaw_scatter_corr = combineScatterSimulation( sinoRaw, primary_bowtie, scatter_bowtie, 10, 1 );

figure; imshow( log10( squeeze( sinoRaw(end/2,:,1:10:end) ) ), [0 6]); title 'Primary photons';
figure; imshow( log10( squeeze( sinoPrime(end/2,:,:)) ), [0 6] ); title 'Primary photons from MC';
figure; imshow( log10( squeeze( sinoScatter(end/2,:,:)) ), [0 6] ); title 'Scattered photons from MC';


%% Process the raw data

geom_orignal = geom;

% Flat filed normalization and log.  
sinoAtt = processCTRawData( sinoRaw, spectrum, tubeCurrentProfile);
% Truncation artifact correction
[ sinoAtt, geom ] = truncationCorrectionWaterCylinderFitting( sinoAtt, geom_orignal, 0.2, 4, 64 );

% same here
sinoAtt_scatter = processCTRawData( sinoRaw_scatter, spectrum, tubeCurrentProfile);
sinoAtt_scatter = truncationCorrectionWaterCylinderFitting( sinoAtt_scatter, geom_orignal, 0.2, 4, 64 );

% same here
sinoAtt_scatter_corr = processCTRawData( sinoRaw_scatter_corr, spectrum, tubeCurrentProfile);
sinoAtt_scatter_corr = truncationCorrectionWaterCylinderFitting( sinoAtt_scatter_corr, geom_orignal, 0.2, 4, 64 );

%% Reconstruction for no scattering case

% setting for reconstruction
kernel = 'hamming';
segmentLength = 'short';
weightingMethod = 'wfdk';

% call reconstruction function
img_wfdk_noscatter = reconHelical( sinoAtt, geom, kernel, weightingMethod, segmentLength, false);

% display
showPhaserRecons( img_wfdk_noscatter, map.windowAtt );


%% Reconstruction for scattering no correction case

% call reconstruction function
img_wfdk_scatter = reconHelical( sinoAtt_scatter, geom, kernel, weightingMethod, segmentLength, false);

% display
showPhaserRecons( img_wfdk_scatter, map.windowAtt );


%% Reconstruction for scattering with perfect correction case

% call reconstruction function
img_wfdk_scatter_corr = reconHelical( sinoAtt_scatter_corr, geom, kernel, weightingMethod, segmentLength, false);

% display
showPhaserRecons( img_wfdk_scatter_corr, map.windowAtt );


