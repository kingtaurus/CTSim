%% test_sys_helical_scatter.m
% This test script is mainly for testing the helical system geometry, 
% helical reconstruction, and scatter simulation.
%
% Meng Wu at Stanford University
% 2013 - 2014

%% Load simulation parameters and datas

load 'temp.mat';

% load phantom
%[phan, map] = loadXCATPhantom(p);
[phan, map ] = loadMaterialsDensityPhantom( p );

%turns = 6;
%pitch = 31/32;
turns = 9;
%pitch = 0.25;
pitch = 0.75;

phan.offset(1) = 30;
phan.offset(2) = 0;

sinosDirKeV = [sinosDirKeV sprintf( '%i-turns-%i-pitch', turns, round( pitch * 100 ) ), '/'];

% Geometry info
geom = loadProjectionGeometryHelicalCT( p, turns, pitch );

% Specturm info
%  spectrum = loadSpectraCT(p, geom, 0.0625*1e6);  % For 1 by 1 detector
  spectrum = loadSpectraCT(p, geom, 0.25*1e6);  % For 2 by 2 detector
% spectrum = loadSpectraCT(p, geom, 1e6);  % For 4 by 4 detector

% Compute ground truth
[imgGtAtt, imgGtHu ] = computeGroundTruth(phan, spectrum);


%% Simulate the raw data

% sinoRaw is the data set without any scatter
[sinoRaw, tubeCurrentProfile] = simulateCTRawData( phan, geom, spectrum, sinosDirKeV );

%% Load MC scatter simulation results

% load( 'E:\MATLAB\CTData\ScatterSimulation\helical\scatter_helical_bowtie_pitch1_32rows.mat' );
% load( 'E:\MATLAB\CTData\scatter_simulation_pelvis\Helical_Pelvis_Final.mat');
load( 'E:\MATLAB\CTData\scatter_simulation_pelvis\Simulation2_Final.mat');

% all_primary = all_primary(:, end/3+1:end*2/3, :);
% all_scatter = all_scatter(:, end/3+1:end*2/3, :);

% all_primary = all_primary(:, 6:end-5, :);
% all_scatter = all_scatter(:, 6:end-5, :);

% sinoRaw_scatter is the data set after adding the simulated scatter, no
% algorithm correction at all (anti-scatter grid is used)

%[sinoRaw_scatter, sinoPrime, sinoScatter] = combineScatterSimulation( sinoRaw, all_primary, all_scatter, 10  );
% sinoRaw_scatter_corr is the data set after adding the simulated scatter,
% no algorithm correction at all (anti-scatter grid is used) 

%sinoRaw_scatter_corr = combineScatterSimulation( sinoRaw, all_primary, all_scatter, 10, 1 );
[sinoRaw_scatter_corr, sinoPrime, sinoScatter_remain] = combineScatterSimulationNew2( sinoRaw, all_primary, all_scatter, 10, 1 );

figure; imshow( log10( squeeze( sinoRaw(end/2,:,1:10:end) ) ), [0 6]); title 'Primary photons';
figure; imshow( log10( squeeze( sinoPrime(end/2,:,:)) ), [0 6] ); title 'Primary photons from MC';
% figure; imshow( log10( squeeze( sinoScatter(end/2,:,:)) ), [0 6] ); title 'Scattered photons from MC';
figure; imdisp( ( sinoScatter_remain(:,:,1:50:end) ), [-1e2 1e2] ); title 'Remained scattered photons after correction';
%figure; imshow( squeeze( sinoRaw_scatter_corr(:,:,1:10:end)) , [0 6] ); title 'Scatter corrected photons from MC';

%% For 256 by 16 detectors, load MC scatter as follows
load( 'E:\MATLAB\CTData\scatter_simulation_pelvis\Simulation2_Final.mat');
sinoRaw = sinoRaw(17:32, :, :);
all_primary = all_primary(:, 17:32, :);
all_scatter = all_scatter(:, 17:32, :);

[sinoRaw_scatter_corr, sinoPrime, sinoScatter_remain] = combineScatterSimulationNew2( sinoRaw, all_primary, all_scatter, 10, 1 );

figure; imshow( log10( squeeze( sinoRaw(end/2,:,1:10:end) ) ), [0 6]); title 'Primary photons';
figure; imshow( log10( squeeze( sinoPrime(end/2,:,:)) ), [0 6] ); title 'Primary photons from MC';
figure; imdisp( ( sinoScatter_remain(:,:,1:50:end) ), [-1e2 1e2] ); title 'Remained scattered photons after correction';

geom.pitch = 0.7500;
geom.detSize = [256 16];

spectrum = loadSpectraCT(p, geom, 1e6);  % For 4 by 4 detector


%% Process the raw data

geom_orignal = geom;

% Flat filed normalization and log.  
sinoAtt = processCTRawData( sinoRaw, spectrum, tubeCurrentProfile);
% Truncation artifact correction
%[ sinoAtt, geom ] = truncationCorrectionWaterCylinderFitting( sinoAtt, geom_orignal, 0.2, 4, 64 );
[ sinoAtt, geom ] = truncationCorrectionWaterCylinderFitting( sinoAtt, geom_orignal, 0.2, 4, 64*4 );

% same here
% sinoAtt_scatter = processCTRawData( sinoRaw_scatter, spectrum, tubeCurrentProfile);
% sinoAtt_scatter = truncationCorrectionWaterCylinderFitting( sinoAtt_scatter, geom_orignal, 0.2, 4, 64 );

% same here
%%sinoAtt_scatter_corr = processCTRawData( sinoRaw_scatter_corr, spectrum, tubeCurrentProfile);
%%sinoAtt_scatter_corr = truncationCorrectionWaterCylinderFitting( sinoAtt_scatter_corr, geom_orignal, 0.2, 4, 64 );

%% Reconstruction for no scattering case

% setting for reconstruction
kernel = 'hamming';
segmentLength = 'short';
weightingMethod = 'wfdk';

% call reconstruction function
img_wfdk_noscatter = reconHelical( sinoAtt, geom, kernel, weightingMethod, segmentLength, false);

% display
showPhaserRecons( img_wfdk_noscatter, map.windowAtt );

imgHu_noscatter = convertMonoAttToHu( img_wfdk_noscatter, spectrum) ;
figure; imdisp( imgHu_noscatter(:,:,end/2), map.windowHu ); 


%% Reconstruction for scattering no correction case

% call reconstruction function
img_wfdk_scatter = reconHelical( sinoAtt_scatter, geom, kernel, weightingMethod, segmentLength, false);

% display
showPhaserRecons( img_wfdk_scatter, map.windowAtt );

imgHu_scatter = convertMonoAttToHu( img_wfdk_noscatter, spectrum) ;
figure; imdisp( imgHu_scatter(:,:,end/2)', map.windowHu ); 


%% Reconstruction for scattering with scatter correction case

% setting for reconstruction
kernel = 'hamming';
segmentLength = 'short';
weightingMethod = 'wfdk';

% call reconstruction function
img_wfdk_scatter_corr = reconHelical( sinoAtt_scatter_corr, geom, kernel, weightingMethod, segmentLength, false);

% display
showPhaserRecons( img_wfdk_scatter_corr, map.windowAtt );

figure; imdisp( img_wfdk_scatter_corr(:,:,round(end/2)), map.windowAtt );

imgHu_corr = convertMonoAttToHu( img_wfdk_scatter_corr, spectrum) ;
figure; imdisp( imgHu_corr(:,:,end/2), map.windowHu ); 



