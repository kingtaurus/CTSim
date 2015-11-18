%% test_sys_helical_scatter_corr.m
% This test script is mainly for testing scatter correction algorithm.
%
% Meng Wu at Stanford University
% 2013 - 2014

%% Load simulation parameters and datas

load 'temp.mat';

% load phantom
[phan, map] = loadXCATPhantom(p);


phan.offset(3) = -93;

turns = 8;
pitch = 0.5;
sinosDirKeV = [sinosDirKeV sprintf( '%i-turns-%i-pitch', turns, round( pitch * 100 ) ), '/'];

% Geometry info
geom = loadProjectionGeometryHelicalCT( p, turns, pitch );

% Specturm info
spectrum = loadSpectraCT(p, geom, 1e6);

% Compute ground truth
[imgGtAtt, imgGtHu ] = computeGroundTruth(phan, spectrum);


% Simulate the raw data

% sinoRaw is the data set without any scatter
[sinoRaw, tubeCurrentProfile] = simulateCTRawData( phan, geom, spectrum, sinosDirKeV, 0 );

%% Load MC scatter simulation results

load( 'E:\MATLAB\CTData\ScatterSimulation\helical\scatter_helical_bowtie_pitch.5_64rows.mat' );

primary_bowtie = rotateSinogram( primary_bowtie, 0, 1, 0 );
primary_bowtie = reverseSingoram( primary_bowtie );


figure; imdisp( log10(sinoRaw(:,:, 1:72:end)) ); title 'Primary photons from MC';
figure; imdisp( log10(primary_bowtie(33:64,:,1:72:end))); title 'Primary photons from MC';


return;

%%
figure; imdisp( log10(scatter_bowtie(:,:,250:260)), [-2 2]); title 'Scattered photons from MC';

% sinoRaw_scatter is the data set after adding the simulated scatter, no
% algorithm correction at all (anti-scatter grid is used)
[sinoRaw_scatter, sinoPrime, sinoScatter] = combineScatterSimulation( sinoRaw, primary_bowtie, scatter_bowtie, 10  );



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
img_wfdk_noscatter = rotateSinogram( img_wfdk_noscatter, 0, 1 );
figure; imdisp( img_wfdk_noscatter(:,:,end:-8:1), map.windowAtt);
figure; imdisp( fliplr( squeeze(img_wfdk_noscatter(:,end/2,:))) , map.windowAtt );
figure; imdisp( fliplr( squeeze(img_wfdk_noscatter(end/2,:,:))) , map.windowAtt );

%% Reconstruction for scattering no correction case

% call reconstruction function
img_wfdk_scatter = reconHelical( sinoAtt_scatter, geom, kernel, weightingMethod, segmentLength, false);

% display
img_wfdk_scatter = rotateSinogram( img_wfdk_scatter, 0, 1 );
figure; imdisp( img_wfdk_scatter(:,:,end:-8:1), map.windowAtt);
figure; imdisp( fliplr( squeeze(img_wfdk_scatter(:,end/2,:))) , map.windowAtt );
figure; imdisp( fliplr( squeeze(img_wfdk_scatter(end/2,:,:))) , map.windowAtt );


%% Reconstruction for scattering with perfect correction case

% call reconstruction function
img_wfdk_scatter_corr = reconHelical( sinoAtt_scatter_corr, geom, kernel, weightingMethod, segmentLength, false);

% display
img_wfdk_scatter_corr = rotateSinogram( img_wfdk_scatter_corr, 0, 1 );
figure; imdisp( img_wfdk_scatter_corr(:,:,end:-8:1), map.windowAtt);
figure; imdisp( fliplr( squeeze(img_wfdk_scatter_corr(:,end/2,:))) , map.windowAtt );
figure; imdisp( fliplr( squeeze(img_wfdk_scatter_corr(end/2,:,:))) , map.windowAtt );

return;
%% In case you get raw scatter simulation data

load( 'E:\MATLAB\CTData\ScatterSimulation\helical\scatter_helical_bowtie_pitch.5_64rows_all.mat' );

primary_bowtie = all_primary_1_2mmAl;
primary_bowtie = primary_bowtie + all_primary_2_3mmAl;
primary_bowtie = primary_bowtie + all_primary_3_4mmAl;
primary_bowtie = primary_bowtie + all_primary_4_5mmAl;
primary_bowtie = primary_bowtie + all_primary_5_6mmAl;
primary_bowtie = primary_bowtie + all_primary_6_7mmAl;
primary_bowtie = primary_bowtie + all_primary_7_8mmAl;
primary_bowtie = primary_bowtie + all_primary_8_9mmAl;
primary_bowtie = primary_bowtie + all_primary_9_10mmAl;
primary_bowtie = primary_bowtie + all_primary_10mmAl;


scatter_bowtie = all_scatter_1_2mmAl;
scatter_bowtie = scatter_bowtie + all_scatter_2_3mmAl;
scatter_bowtie = scatter_bowtie + all_scatter_3_4mmAl;
scatter_bowtie = scatter_bowtie + all_scatter_4_5mmAl;
scatter_bowtie = scatter_bowtie + all_scatter_5_6mmAl;
scatter_bowtie = scatter_bowtie + all_scatter_6_7mmAl;
scatter_bowtie = scatter_bowtie + all_scatter_7_8mmAl;
scatter_bowtie = scatter_bowtie + all_scatter_8_9mmAl;
scatter_bowtie = scatter_bowtie + all_scatter_9_10mmAl;
scatter_bowtie = scatter_bowtie + all_scatter_10mmAl;

figure; imdisp( log10( primary_bowtie(:,:,250) ), [0 1])
figure; imdisp( log10( scatter_bowtie(:,:,250) ), [0 1])

