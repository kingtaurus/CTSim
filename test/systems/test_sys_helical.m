%% test_sys_helical.m
% This test script is mainly for testing the helical system geometry, 
% helical reconstruction, and cone-beam artifacts
%
% Meng Wu at Stanford University
% 2013 - 2014

%% Load simulation parameters and datas

load 'temp.mat';

% load phantom
[phan, map] = loadXCATPhantom(p);

turns = 8;
% feel free to change the helical pitch
pitch = 0.5; 

sinosDir = [sinosDir sprintf( '%i-turns-%i-pitch', turns, round( pitch * 100 ) ), '/'];

% Geometry info
geom = loadProjectionGeometryHelicalCT( p, turns, pitch );

% Specturm info
spectrum = loadSpectraCT(p, geom, 2e5);

% Compute ground truth
[imgGtAtt, imgGtHu ] = computeGroundTruth(phan, spectrum);

%% Simulate the raw data

[sinoRaw, tubeCurrentProfile] = simulateCTRawData( phan, geom, spectrum, sinosDirKeV );

%% Process the raw data. 
% Flat filed normalization and log.  

sinoAtt = processCTRawData( sinoRaw, spectrum, tubeCurrentProfile);

%% Truncation artifact correction

geom_orignal = geom;
[ sinoAtt_trun_corr, geom_trun_corr ] = truncationCorrectionWaterCylinderFitting( sinoAtt, geom_orignal, 0.2, 4, 64 );

%% Weighted FDK reconstruction without trucation correction

% setting for reconstruction
kernel = 'hamming';
segmentLength = 'short';
weightingMethod = 'wfdk';

% call reconstruction function
img_wfdk = reconHelical( sinoAtt, geom_orignal, kernel, weightingMethod, segmentLength, false);

% display
showPhaserRecons( img_wfdk, map.windowAtt );

%% Weighted FDK reconstruction with trucation correction

% call reconstruction function
img_wfdk_trun_corr = reconHelical( sinoAtt_trun_corr, geom_trun_corr, kernel, weightingMethod, segmentLength, false);

% display
showPhaserRecons( img_wfdk_trun_corr, map.windowAtt );


%% View-weighted FDK reconstruction with trucation correction

% change weighting method to 'vwfdk'
img_vwfdk_trun_corr = reconHelical( sinoAtt_trun_corr, geom_trun_corr, kernel, 'vwfdk', segmentLength, false);

% display
showPhaserRecons( img_vwfdk_trun_corr, map.windowAtt );


%% Two-pass cone beam artifact correction version 1

% initalize baseline recon method
reconFunc = @(sino)reconHelical( sino, geom_trun_corr, kernel, weightingMethod, segmentLength);

% call reconstruction function
img_wfdk_2pass_1 = reconTwoPassConebeamArtifactCorrectionBone( sinoAtt_trun_corr, geom_trun_corr, reconFunc, 0.2 );

% display
showPhaserRecons( img_wfdk_2pass_1, map.windowAtt );

%% Two-pass cone beam artifact correction version 2

% call reconstruction function
img_wfdk_2pass_2 = reconTwoPassConebeamArtifactCorrectionTissueBone( sinoAtt_trun_corr, geom_trun_corr, reconFunc, 0.2 );

% display
showPhaserRecons( img_wfdk_2pass_2, map.windowAtt );

