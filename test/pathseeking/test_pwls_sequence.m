%load 'temp.mat';

% Load simulation parameters and datas
[phan, map] = loadXCATPhantom(p);

[ geom ] = loadProjectionGeometryCT( p );

spectrum = loadSpectraCT(p, geom, 1e6);

% compute sinogram

[sinoRaw, tubeCurrentProfile] = simulateCTRawData( phan, geom, spectrum, sinosDir, 1, 1, 2 );

sinoAtt = processCTRawData( sinoRaw, spectrum, tubeCurrentProfile);

%Compute ground trut
[imgGtAtt, imgGtHu ] = computeGroundTruth(phan, spectrum);

weights = computeWeightsPwls( sinoRaw, 0, spectrum.electronicNoise );

%%
nitn = 50;
numos = 8;
beta1 = 1e3;
beta2 = 1e5;
delta = 1e-3;

% Nesterov + OS + SQS
img_fbp = reconFBP( sinoAtt, geom, 'ram-lak');
impop( img_fbp );


[img_pwls1] = reconPwlsSeNesterovSqs( sinoAtt, weights, geom, beta1, 'huber', nitn, delta, numos );
impop( img_pwls1 );

[img_pwls2] = reconPwlsSeNesterovSqs( sinoAtt, weights, geom, beta2, 'huber', nitn, delta, numos );
impop( img_pwls2 );

img_true = reconPwlsSequencesNesterov( sinoAtt, weights, geom, [beta1 beta2], 'huber', delta, numos, 10 );
figure; imdisp( img_true(:,:,1:10), [0.16 0.26] );

