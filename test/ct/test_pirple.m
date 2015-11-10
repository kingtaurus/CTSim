%% test_sparse_view.m
% Compare several reconstruction methods for sparse view CT
% Meng Wu at University of Erlangen-Nuremburg
% 2014.10

%% Simulate normal CT data
clear;
close all;
load 'temp.mat';

% Load simulation parameters and datas
[phan, map] = loadXCATPhantom(p);

geom = loadProjectionGeometryShortScan( p, 223, 90 );

% geom = loadProjectionGeometryCT( p );

spectrum = loadSpectraCT(p, geom, 1e5);

% compute sinogram
[sinoRaw, tubeCurrentProfile] = simulateCTRawData( phan, geom, spectrum, sinosDir, 1, 1, 2 );

sinoAtt = processCTRawData( sinoRaw, spectrum, tubeCurrentProfile);

weights = computeWeightsPwls( sinoRaw, 0, spectrum.electronicNoise );

%% Reconstruction without view downsampling 

% Filter backprojection with hamming window

img_fbp = reconFBP( sinoAtt, geom, 'hamming', 1, 1);

figure; imshow( img_fbp(:,:,ceil(end/2)), map.windowAtt ); 
title 'FBP reconstruction hamming window (full views)';
evaluateSoftTissue( img_fbp, spectrum, map );

%% Downsample projection data to mimic sparse view case 
% We first consider downsample by 4 which is a moderate sparse view case,
% hard enough for interpolation and FBP reconstruction. The TV iterative
% reconstruction can deal it well.

rate = 36;

[geomSparse, sinoSparse, weightsSparse ]  = convertSinogram2SparseViewScan( geom, sinoAtt, rate, weights );

img_sparse_fbp = reconFBP( sinoSparse, geomSparse, 'hamming', 1, 1 );

figure; imshow( img_sparse_fbp(:,:,ceil(end/2)), map.windowAtt ); 
title 'FBP reconstruction (sparse views)';
evaluateSoftTissue( img_sparse_fbp, spectrum, map );

%%

nitn = 20;
numos = 1;
beta = 200;
delta = 1e-3;

img0 =  reconFBP( sinoSparse, geomSparse, 'hamming', 1, 0.3 );

img_sparse_tv = reconPwlsLALMOs14( sinoSparse, weightsSparse, geomSparse, beta, 'isotv', nitn, delta, numos, img0  );

figure; imshow( img_sparse_tv(:,:,ceil(end/2)), map.windowAtt );
title( ['PWLS with TV reconstruction beta = ' num2str(beta)] );
evaluateSoftTissue( img_sparse_tv, spectrum, map );

