load 'temp.mat';

% Load simulation parameters and datas
[phan, map] = loadXCATPhantom(p);

[ geom ] = loadProjectionGeometryCT( p );

spectrum = loadSpectraCT(p, geom, 1e5);

% compute sinogram
sinoPhotonCounts= simulatePhotonCountingData( phan, geom, spectrum, sinosDirKeV );
sinoAtt = computePhotonCountingSinogramAttunation( sinoPhotonCounts, spectrum );

% corrections
sinoAtt = beamHardeningWarterCorrection(sinoAtt, spectrum);

%Compute ground trut
[imgGtAtt, imgGtHu ] = computeGroundTruth(phan, spectrum);

weights = computeWeightsPwls( sinoPhotonCounts, 0, spectrum.backgroundEvents );
nitn = 20;
beta = 1e4;
delta = 1e-3;

%load 'pwls_conv_2.mat';

img_converged = reconPwlsSeNesterovSqs( sinoAtt, weights, geom, beta, 'huber', 500, delta, 4);
impop( img_converged );

%% FBP

img_fbp = reconFBP( sinoAtt, geom, 'hamming', 1, 1);
impop( img_fbp );

%% OS + SQS

[img_sqs, phis_sqs, rmsds_sqs ] = reconPwlsSeSqs( sinoAtt, weights, geom, beta, 'huber', nitn, delta, 20, img_fbp, img_converged);
impop( img_sqs );

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



