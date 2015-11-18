%% Test new geometry
%  The detectors and sources are on the same cylinder 
%  He Yang, Feb 2015

%% Load parameters and data

load 'temp.mat';

% Load material and density
[phan, map ] = loadMaterialsDensityPhantom( p );

% Load system geometry
[ geom ] = loadProjectionGeometryCT( p );

% Load spectrum
%spectrum = loadSpectraCT(p, geom, 0.25*1e6);
spectrum = loadSpectraCT(p, geom, 1e6);

%% Simulate the raw data

nzSource = 3;
[sinoRaw, tubeCurrentProfile] = simulateCTRawData2( phan, geom, spectrum, sinosDirKeV, nzSource); % New Geometry
%[sinoRaw, tubeCurrentProfile] = simulateCTRawData( phan, geom, spectrum, sinosDirKeV); % Old Geometry

% if nzSource ~=1
%     tubeCurrentProfile = nzSource*tubeCurrentProfile;
% end

sinoAtt = processCTRawData( sinoRaw, spectrum, tubeCurrentProfile);

% compute sinogram
%sinoPhotonCounts= simulatePhotonCountingData( phan, geom, spectrum, sinosDirKeV );
%sinoAtt = computePhotonCountingSinogramAttunation( sinoPhotonCounts, spectrum );

% corrections
%sinoAtt = beamHardeningWarterCorrection(sinoAtt, spectrum);

%Compute ground trut
%[imgGtAtt, imgGtHu ] = computeGroundTruth(phan, spectrum);

weights = computeWeightsPwls( sinoRaw, 0, spectrum.electronicNoise );
%nitn = 20;
%beta = 1e4;
beta = 5e4;
delta = 1e-3;

%%
%load 'pwls_conv_2.mat';

img_converged = reconPwlsSeNesterovSqs( sinoAtt, weights, geom, beta, 'huber', 500, delta, 4);
impop( img_converged );

%% FBP
% 
% img_fbp = reconFBP( sinoAtt, geom, 'hamming', 1, 1);
% figure; imagesc( img_fbp(:,:,round(end/2)), map.windowAtt);


%% OS + SQS
% New geometry (changed loadPojectors3.m to load new forward and backward projectors)

nitn = 40;
beta = 1e5;
%img0 = img_sqs;
img0 = zeros(geom.reconSize);

% [img_sqs, phis_sqs, rmsds_sqs ] = reconPwlsSeSqs3( sinoAtt, weights, geom, beta, 'huber', nitn, delta, 20, img_fbp, img_converged);
%impop( img_sqs );
%nzSource = 1;

[img_sqs, phis_sqs, rmsds_sqs ] = reconPwlsSeSqs3( sinoAtt, weights, geom, beta, 'huber', nitn, delta, 20, img0, nzSource);
%figure; imagesc( img_sqs(:,:,round(end/2)), map.windowAtt);

imgHu_sqs = convertMonoAttToHu( img_sqs, spectrum); 
%figure; imdisp( imgHu_sqs(:,:,round(end/2)), map.windowHu);
figure; imdisp( imgHu_sqs(:,:,round(end/2)), [-1000 1000]);
%figure; imdisp( (imgHu_sqs(:,:,round(end/2)) + imgHu_sqs(:,:,round(end/2)+1))/2, [-1000 1000]);
%figure; imdisp( (imgHu_sqs(:,:,round(end/2)) + imgHu_sqs(:,:,round(end/2)+1))/2 - imgGtHu, [-1000 1000]);
figure; imdisp( imgHu_sqs(:,end/2,:), [-1000 1000]);


% Compute ground truth
spectrum1 = spectrum;
spectrum1.energyAverage = 72;
[imgGtAtt, imgGtHu ] = computeGroundTruth(phan, spectrum1);

figure; imdisp( imgGtHu, [-1000 1000]);
figure; imdisp( imgHu_sqs(:,:,round(end/2)) - imgGtHu, [-1000 1000]);


%% OS + NU-SQS

[img_nusqs, phis_nusqs, rmsds_nusqs ] = reconPwlsSeNusqs( sinoAtt, weights, geom, beta, 'huber', nitn, delta, 20, img_fbp, img_converged);
impop( img_nusqs );

%% Nesterov + OS + SQS

[img_nv_sqs, phis_nv_sqs, rmsds_nv_sqs ] = reconPwlsSeNesterovSqs( sinoAtt, weights, geom, beta, 'huber', nitn, delta, 10, img_fbp, img_converged );
impop( img_nv_sqs );

%% Nesterov + OS + NU-SQS

[img_nv_nusqs, phis_nv_nusqs, rmsds_nv_nusqs ] = reconPwlsSeNesterovNusqs( sinoAtt, weights, geom, beta, 'huber', nitn, delta, 10, img_fbp, img_converged);
impop( img_nv_nusqs );

%% ADMM 

[img_admm, phis_admm, rmsds_admm ] = reconPwlsADMM( sinoAtt, weights, geom, beta, 'huber', nitn, delta, 1, img_fbp, img_converged);
impop( img_admm );

%% LALM

[img_lalm, phis_lalm, rmsds_lalm ] = reconPwlsLALM( sinoAtt, weights, geom, beta, 'huber', nitn, delta, 1, img_fbp, img_converged);
impop( img_lalm );

%% LALM + OS

[img_lalm_os, phis_lalm_os, rmsds_lalm_os ] = reconPwlsLALMOs( sinoAtt, weights, geom, beta, 'huber', nitn, delta, 10, img_fbp,  img_converged);
impop( img_lalm_os );


%% LALM + OS + 14

[img_lalm_os14, phis_lalm_os14, rmsds_lalm_os14 ] = reconPwlsLALMOs14( sinoAtt, weights, geom, beta, 'huber', nitn, delta, 10, img_fbp, img_converged);
impop( img_lalm_os14 );



