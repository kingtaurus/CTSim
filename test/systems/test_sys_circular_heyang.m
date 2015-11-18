%% test_sys_circular_heyang.m
% This test script is for testing the circular system geometry, 
% circular reconstruction, and scatter simulation 
% (based on test_sys_helical_scatter.m)
%
% Last updated by He Yang
% Nov 2014

%% Load simulation parameters and data

load 'temp.mat';

% Load phantom
%[phan, map] = loadXCATPhantom(p);
[phan, map ] = loadMaterialsDensityPhantom( p );
% Load distorted pelvis phantom
%[phan, map] = loadXCATPhantomDistorted(p);

%sinosDirKeV = [sinosDirKeV sprintf( '%i-turns-%i-pitch', turns, round( pitch * 100 ) ), '/'];

% phan.offset(1) = 25;
% phan.offset(2) = 0;

% Geometry info
geom = loadProjectionGeometryCT( p );

% Spectrum info
spectrum = loadSpectraCT(p, geom, 1e6);  % For 4 by 4 detector
% spectrum = loadSpectraCT(p, geom, 0.25*1e6); % For 2 by 2 detector
% spectrum = loadSpectraCT(p, geom, 0.0625*1e6); % For 1 by 1 detector
%sinoAtt = simulateAttenuationDataPhantom( phan, geom, spectrum, sinosDirKeV);

% Compute ground truth
[imgGtAtt, imgGtHu ] = computeGroundTruth(phan, spectrum);


%% Simulate the raw data

% sinoRaw is the data set without any scatter
[sinoRaw, tubeCurrentProfile] = simulateCTRawData( phan, geom, spectrum, sinosDirKeV);

%% Load MC scatter simulation results

% load( 'E:\MATLAB\CTData\scatter_simulation_pelvis\PelvisData_Final.mat');
% load( 'E:\MATLAB\CTData\scatter_simulation_pelvis\Pelvis3_Cone_Scatter.mat');
% load( 'E:\MATLAB\CTData\scatter_simulation_pelvis\Pelvis5_Cone.mat');
%load( 'E:\MATLAB\CTData\scatter_simulation_pelvis\Pelvis5_Cone_corrected.mat');
load( 'E:\MATLAB\CTData\scatter_simulation_pelvis\Simulation5_Final.mat');
% bin 4*4 into 2*2
% all_primary = phaserLargeDetectorBinning2(all_primary);
% all_scatter = phaserLargeDetectorBinning2(all_scatter);

% sinoRaw_scatter is the data set after adding the simulated scatter
% no algorithm correction at all (anti-scatter grid is used)

%[sinoRaw_scatter, sinoPrime, sinoScatter] = combineScatterSimulation( sinoRaw, primary_final, scatter_final, 10);
%[sinoRaw_scatter, sinoPrime, sinoScatter] = combineScatterSimulation( sinoRaw, all_primary, all_scatter, 10);
%[sinoRaw_scatter, sinoPrime, sinoScatter] = combineScatterSimulation( sinoRaw, all_primary, all_scatter, 10);
% sinoRaw_scatter_corr is the data set after adding the simulated scatter,
% no algorithm correction at all (anti-scatter grid is used) 

%sinoRaw_scatter_corr = combineScatterSimulation( sinoRaw, primary_final, scatter_final, 10, 1 );
%sinoRaw_scatter_corr = combineScatterSimulation( sinoRaw, all_primary, all_scatter, 10, 1 );
%[sinoRaw_scatter_corr, sinoPrime, sinoScatter] = combineScatterSimulationNew( sinoRaw, all_primary, all_scatter, 10,1); % With simple scatter correction
[ sinoRaw_scatter_corr, sinoPrime, sinoScatter_remain ] = combineScatterSimulationNew( sinoRaw, all_primary, all_scatter, 10,1); % With simple scatter correction
%sinoRaw_scatter_corr = combineScatterSimulation( sinoRaw, primary, scatter, 90, 1 );

%figure; imshow( log10( squeeze( sinoRaw(end/2,:,1:10:end) ) ), [0 6]); title 'Primary photons';
figure; imdisp( ( squeeze( sinoPrime(:,:,:)) ), [0 1e4] ); title 'Primary photons from MC';
figure; imdisp( ( sinoScatter_remain(:,:,:) ), [-1e2 1e2] ); title 'Scattered photons from MC';
figure; imdisp( ( squeeze( sinoRaw_scatter_corr(:,:,1:10:end) ) ), [0 1e4] ); title 'Scatter corrected primary';

%% Process the raw data

geom_orignal = geom;

% Takes uniform detector and bins outer edge as 4*4 to account for larger pixels near the outside 
%sinoRaw = phaserLargeDetectorBinning( sinoRaw );
sinoRaw = phaserLargeDetectorBinning2( sinoRaw ); % bins 4*4 pixels into 2*2 pixels

% Flat filed normalization and log.  
sinoAtt = processCTRawData( sinoRaw, spectrum, tubeCurrentProfile);
% Truncation artifact correction
[ sinoAtt, geom ] = truncationCorrectionWaterCylinderFitting( sinoAtt, geom_orignal, 0.2, 4, 64 );

% same here
% sinoAtt_scatter = processCTRawData( sinoRaw_scatter, spectrum, tubeCurrentProfile);
% sinoAtt_scatter = truncationCorrectionWaterCylinderFitting( sinoAtt_scatter, geom_orignal, 0.2, 4, 64 );

% same here
sinoAtt_scatter_corr = processCTRawData( sinoRaw_scatter_corr, spectrum, tubeCurrentProfile);
sinoAtt_scatter_corr = truncationCorrectionWaterCylinderFitting( sinoAtt_scatter_corr, geom_orignal, 0.2, 4, 64 );


%% Padding the z-direction
 [nu, nv, nview] = size(sinoAtt_scatter_corr);
sinoAtt_scatter_corr1 =  zeros( nu+4*2, nv , nview, 'single' ); % Extend in z-direction
sinoAtt_scatter_corr1(5:52, 1:nv, 1:nview) = sinoAtt_scatter_corr;
sinoAtt_scatter_corr1(1:4,:,:) = repmat(sinoAtt_scatter_corr(1,:,:), 4, 1);
sinoAtt_scatter_corr1(end-3:end,:,:) = repmat(sinoAtt_scatter_corr(end,:,:), 4, 1);

geom.detSize(2) = geom.detSize(2) + 4*2;


%% Reconstruction for no scattering case

% setting for reconstruction
%kernel = 'hamming';
%segmentLength = 'short';
%weightingMethod = 'wfdk';
map.windowHu = [-1000 1000];

% call reconstruction function
imgAttFBP_noscatter = reconFBP( sinoAtt, geom, 'hamming', 1, 1);
%img_wfdk_noscatter = reconHelical( sinoAtt, geom, kernel, weightingMethod, segmentLength, false);

% display
% img_wfdk_noscatter = rotateSinogram( img_wfdk_noscatter, 0, 1 );
% figure; imdisp( img_wfdk_noscatter(:,:,end:-8:1), map.windowAtt);
% figure; imdisp( fliplr( squeeze(img_wfdk_noscatter(:,end/2,:))) , map.windowAtt );
% figure; imdisp( fliplr( squeeze(img_wfdk_noscatter(end/2,:,:))) , map.windowAtt );

imgAttFBP_noscatter = rotateSinogram( imgAttFBP_noscatter, 0, 1 );
imgFBP_noscatter = convertMonoAttToHu( imgAttFBP_noscatter, spectrum) ;
figure; imdisp( imgFBP_noscatter(:,:,end/2)', map.windowHu ); 
%figure; imdisp( fliplr( squeeze(imgFBP_noscatter(end/2,:,:))), map.windowHu ); 
%figure; imdisp( imgFBP_noscatter(:,:,:), map.windowHu ); 
figure; imdisp( fliplr( squeeze(imgFBP_noscatter(:,end/2,:) ) ), map.windowHu );
%figure; imdisp( fliplr( squeeze(imgFBP_noscatter(end/2,:,:) ) ), map.windowHu );
figure; imdisp( imgFBP_noscatter(:,:,:) , map.windowHu );


%% Reconstruction for scattering no correction case

% call reconstruction function
imgAttFBP_scatter = reconFBP( sinoAtt_scatter, geom, 'hamming', 1, 1);

% display
imgAttFBP_scatter = rotateSinogram( imgAttFBP_scatter, 0, 1 );
imgFBP_scatter = convertMonoAttToHu( imgAttFBP_scatter, spectrum) ;
figure; imdisp( imgFBP_scatter(:,:,end/2)', map.windowHu ); 
%figure; imdisp( fliplr( squeeze(imgFBP_scatter(end/2,:,:))), map.windowHu ); 
%figure; imdisp( imgFBP_scatter(:,:,:), map.windowHu ); 
figure; imdisp( fliplr( squeeze(imgFBP_scatter(:,end/2,:) ) ), map.windowHu );
%figure; imdisp( fliplr( squeeze(imgFBP_scatter(end/2,:,:) ) ), map.windowHu );
figure; imdisp( imgFBP_scatter(:,:,:), map.windowHu );

%% Reconstruction for scattering with scatter correction

map.windowHu = [-1000 1000];

% call reconstruction function
imgAttFBP_scatter_corr = reconFBP( sinoAtt_scatter_corr1, geom, 'hamming',1,1);

% display
imgAttFBP_scatter_corr = rotateSinogram( imgAttFBP_scatter_corr, 0, 1 );
imgFBP_scatter_corr = convertMonoAttToHu( imgAttFBP_scatter_corr, spectrum) ;
figure; imdisp( imgFBP_scatter_corr(:,:,end/2)', map.windowHu ); 
%figure; imdisp( fliplr( squeeze(imgFBP_scatter_corr(end/2,:,:))), map.windowHu ); 
%figure; imdisp( imgFBP_scatter_corr(:,:,:), map.windowHu ); 
figure; imdisp( fliplr( squeeze(imgFBP_scatter_corr(:,end/2,:) ) ), map.windowHu );
%figure; imdisp( fliplr( squeeze(imgFBP_scatter_corr(end/2,:,:) ) ), map.windowHu );
figure; imdisp( imgFBP_scatter_corr(:,:,:), map.windowHu ); 



%%
%figure; imdisp( imgAttFBP_scatter_corr(:,:,end:-8:1), map.windowAtt);
%figure; imdisp( fliplr( squeeze(imgAttFBP_scatter_corr(:,end/2,:))) , map.windowAtt );
%figure; imdisp( fliplr( squeeze(imgAttFBP_scatter_corr(end/2,:,:))) , map.windowAtt );

%%
% sinoRaw = phaserLargeDetectorBinning( sinoRaw );
% 
% sinoAtt =  processCTRawData( sinoRaw, spectrum, tubeCurrentProfile);
% 
% geom_orignal = geom;
% [ sinoAtt_trun_corr, geom_trun_corr ] = truncationCorrectionWaterCylinderFitting( sinoAtt, geom_orignal, 0.2, 4, 64 );

% %% full scan
%  map.windowHu  = [-1000 1000];
% 
% 
% imgAttFBP = reconFBP( sinoAtt, geom, 'hamming', 1, 1);
% imgFBP = convertMonoAttToHu( imgAttFBP, spectrum) ;
% figure; imdisp( imgFBP(:,:,end/2)', map.windowHu ); 
% figure; imdisp( fliplr( squeeze(imgFBP(end/2,:,:))), map.windowHu ); 
% figure; imdisp( imgFBP(:,:,:), map.windowHu ); 


%% short scan
offsetAngle = 0;

[geomShort, sinoShort ]  = convertSinogram2ShortScan( geom, sinoAtt, offsetAngle );
imgAttShort = reconFBP( sinoShort, geomShort, 'hamming' );
imgFBPShort = convertMonoAttToHu( imgAttShort, spectrum) ;
figure; imdisp( imgFBPShort(:,:,end/2)', map.windowHu ); 
figure; imdisp( fliplr( squeeze(imgFBPShort(end/2,:,:))), map.windowHu ); 
figure; imdisp( imgFBPShort(:,:,:), map.windowHu ); 

return