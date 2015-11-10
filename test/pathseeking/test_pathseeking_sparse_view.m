%% test_sparse_view.m
% Compare several reconstruction methods for sparse view CT
% Meng Wu at University of Erlangen-Nuremburg
% 2014.10

%% Simulate normal CT data
clear;
close all;
load 'temp.mat';

% Load simulation parameters and datas
[phan, map] = loadMaterialsDensityPhantom(p);

geom = loadProjectionGeometryShortScan( p, 223, 90 );

% geom = loadProjectionGeometryCT( p );

spectrum = loadSpectraCT(p, geom, 1e6);

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

rate = 8;

[geomSparse, sinoSparse, weightsSparse ]  = convertSinogram2SparseViewScan( geom, sinoAtt, rate, weights );

img_sparse_fbp = reconFBP( sinoSparse, geomSparse, 'hamming', 1, 1 );

figure; imshow( img_sparse_fbp(:,:,ceil(end/2)), map.windowAtt ); 
title 'FBP reconstruction (sparse views)';
evaluateSoftTissue( img_sparse_fbp, spectrum, map );


%% Linear interpolation method

[sinoInterp, geomInterp] = upsamplingViews( sinoSparse, geomSparse, rate );

img_sparse_fbp_interp = reconFBP( sinoInterp, geomInterp, 'hamming', 1, 1);

figure; imshow( img_sparse_fbp_interp(:,:,ceil(end/2)), map.windowAtt );
title 'FBP linear interpolation reconstruction';
evaluateSoftTissue( img_sparse_fbp_interp, spectrum, map );

img0 = img_sparse_fbp_interp;


%%
nitn = 100;
numos = 4;
beta1 = 20;
beta2 = 2000;
delta = 1e-3;
noframes = 20;


img_sparse_tv_1 = reconPwlsLALMOs14( sinoSparse, weightsSparse, geomSparse, beta1, 'isotv', nitn, 1, numos, img0  );

figure; imshow( img_sparse_tv_1(:,:,ceil(end/2)), map.windowAtt );
title( ['PWLS with TV reconstruction beta = ' num2str(beta1)] );
evaluateSoftTissue( img_sparse_tv_1, spectrum, map );

img_sparse_tv_2 = reconPwlsLALMOs14( sinoSparse, weightsSparse, geomSparse, beta2, 'isotv', nitn, 1, numos, img0  );

figure; imshow( img_sparse_tv_2(:,:,ceil(end/2)), map.windowAtt );
title( ['PWLS with TV reconstruction beta = ' num2str(beta2)] );
evaluateSoftTissue( img_sparse_tv_2, spectrum, map );

img_seqs = reconPwlsSequencesLALM(  sinoSparse, weightsSparse, geomSparse, [beta1 beta2], 'isotv', delta, numos, nitn, noframes );

figure; imdisp( img_seqs(:,:,:), [0.1 0.3] );



%%

noframes = 40;

dv = 2e-4;
p = 0.2;
os = 3;

img_aps_f_os3_dv1_p2 = reconPwlsApproxPathSeeking( sinoSparse, weightsSparse, geomSparse, 'isotv', delta, os, img_sparse_tv_1, img_sparse_tv_2, dv, p, noframes);

figure; imdisp( img_aps_f_os3_dv1_p2, [0.1 0.3]  );

img_tps1_f_os3_dv1_p2 = reconPwlsTruePathSeekingIndep( sinoSparse, weightsSparse, geomSparse, 'isotv', delta, os, img_sparse_tv_1, img_sparse_tv_2, dv, p, noframes);

figure; imdisp( img_tps1_f_os3_dv1_p2, [0.1 0.3]  );

img_tps2_f_os3_dv1_p2 = reconPwlsTruePathSeekingSelf( sinoSparse, weightsSparse, geomSparse, 'isotv', delta, os, img_sparse_tv_1, img_sparse_tv_2, dv, p, noframes);

figure; imdisp( img_tps2_f_os3_dv1_p2, [0.1 0.3]  );

img_aps_b_os3_dv1_p2 = reconPwlsApproxPathSeeking( sinoSparse, weightsSparse, geomSparse, 'isotv', delta, os, img_sparse_tv_2, img_sparse_tv_1, dv, p, noframes);

figure; imdisp( img_aps_b_os3_dv1_p2, [0.1 0.3]  );

img_tps1_b_os3_dv1_p2 = reconPwlsTruePathSeekingIndep( sinoSparse, weightsSparse, geomSparse, 'isotv', delta, os, img_sparse_tv_2, img_sparse_tv_1, dv, p, noframes);

figure; imdisp( img_tps1_b_os3_dv1_p2, [0.1 0.3]  );

img_tps2_b_os3_dv1_p2 = reconPwlsTruePathSeekingSelf( sinoSparse, weightsSparse, geomSparse, 'isotv', delta, os, img_sparse_tv_2, img_sparse_tv_1, dv, p, noframes);

figure; imdisp( img_tps2_b_os3_dv1_p2, [0.1 0.3]  );

%
dv = 2e-4;
p = 0.3;
os = 3;

img_aps_f_os3_dv1_p3 = reconPwlsApproxPathSeeking( sinoSparse, weightsSparse, geomSparse, 'isotv', delta, os, img_sparse_tv_1, img_sparse_tv_2, dv, p, noframes);

figure; imdisp( img_aps_f_os3_dv1_p3, [0.1 0.3]  );

img_tps1_f_os3_dv1_p3 = reconPwlsTruePathSeekingIndep( sinoSparse, weightsSparse, geomSparse, 'isotv', delta, os, img_sparse_tv_1, img_sparse_tv_2, dv, p, noframes);

figure; imdisp( img_tps1_f_os3_dv1_p3, [0.1 0.3]  );

img_tps2_f_os3_dv1_p3 = reconPwlsTruePathSeekingSelf( sinoSparse, weightsSparse, geomSparse, 'isotv', delta, os, img_sparse_tv_1, img_sparse_tv_2, dv, p, noframes);

figure; imdisp( img_tps2_f_os3_dv1_p3, [0.1 0.3]  );

img_aps_b_os3_dv1_p3 = reconPwlsApproxPathSeeking( sinoSparse, weightsSparse, geomSparse, 'isotv', delta, os, img_sparse_tv_2, img_sparse_tv_1, dv, p, noframes);

figure; imdisp( img_aps_b_os3_dv1_p3, [0.1 0.3]  );

img_tps1_b_os3_dv1_p3 = reconPwlsTruePathSeekingIndep( sinoSparse, weightsSparse, geomSparse, 'isotv', delta, os, img_sparse_tv_2, img_sparse_tv_1, dv, p, noframes);

figure; imdisp( img_tps1_b_os3_dv1_p3, [0.1 0.3]  );

img_tps2_b_os3_dv1_p3 = reconPwlsTruePathSeekingSelf( sinoSparse, weightsSparse, geomSparse, 'isotv', delta, os, img_sparse_tv_2, img_sparse_tv_1, dv, p, noframes);

figure; imdisp( img_tps2_b_os3_dv1_p3, [0.1 0.3]  );

%
dv = 4e-4;
p = 0.2;
os = 3;

img_aps_f_os3_dv2_p2 = reconPwlsApproxPathSeeking( sinoSparse, weightsSparse, geomSparse, 'isotv', delta, os, img_sparse_tv_1, img_sparse_tv_2, dv, p, noframes);

figure; imdisp( img_aps_f_os3_dv2_p2, [0.1 0.3]  );

img_tps1_f_os3_dv2_p2 = reconPwlsTruePathSeekingIndep( sinoSparse, weightsSparse, geomSparse, 'isotv', delta, os, img_sparse_tv_1, img_sparse_tv_2, dv, p, noframes);

figure; imdisp( img_tps1_f_os3_dv2_p2, [0.1 0.3]  );

img_tps2_f_os3_dv2_p2 = reconPwlsTruePathSeekingSelf( sinoSparse, weightsSparse, geomSparse, 'isotv', delta, os, img_sparse_tv_1, img_sparse_tv_2, dv, p, noframes);

figure; imdisp( img_tps2_f_os3_dv2_p2, [0.1 0.3]  );

img_aps_b_os3_dv2_p2 = reconPwlsApproxPathSeeking( sinoSparse, weightsSparse, geomSparse, 'isotv', delta, os, img_sparse_tv_2, img_sparse_tv_1, dv, p, noframes);

figure; imdisp( img_aps_b_os3_dv2_p2, [0.1 0.3]  );

img_tps1_b_os3_dv2_p2 = reconPwlsTruePathSeekingIndep( sinoSparse, weightsSparse, geomSparse, 'isotv', delta, os, img_sparse_tv_2, img_sparse_tv_1, dv, p, noframes);

figure; imdisp( img_tps1_b_os3_dv2_p2, [0.1 0.3]  );

img_tps2_b_os3_dv2_p2 = reconPwlsTruePathSeekingSelf( sinoSparse, weightsSparse, geomSparse, 'isotv', delta, os, img_sparse_tv_2, img_sparse_tv_1, dv, p, noframes);

figure; imdisp( img_tps2_b_os3_dv2_p2, [0.1 0.3]  );

%%

dv = 4e-4;
p = 0.1;
os = 3;

img_aps_f_os3_dv2_p1 = reconPwlsApproxPathSeeking( sinoSparse, weightsSparse, geomSparse, 'isotv', delta, os, img_sparse_tv_1, img_sparse_tv_2, dv, p, noframes);

figure; imdisp( img_aps_f_os3_dv2_p1, [0.15 0.25]  );

img_tps1_f_os3_dv2_p1 = reconPwlsTruePathSeekingIndep( sinoSparse, weightsSparse, geomSparse, 'isotv', delta, os, img_sparse_tv_1, img_sparse_tv_2, dv, p, noframes);

figure; imdisp( img_tps1_f_os3_dv2_p1, [0.1 0.3]  );

img_tps2_f_os3_dv2_p1 = reconPwlsTruePathSeekingSelf( sinoSparse, weightsSparse, geomSparse, 'isotv', delta, os, img_sparse_tv_1, img_sparse_tv_2, dv, p, noframes);

figure; imdisp( img_tps2_f_os3_dv2_p1, [0.15 0.25]  );

img_aps_b_os3_dv2_p1 = reconPwlsApproxPathSeeking( sinoSparse, weightsSparse, geomSparse, 'isotv', delta, os, img_sparse_tv_2, img_sparse_tv_1, dv, p, noframes);

figure; imdisp( img_aps_b_os3_dv2_p1, [0.1 0.3]  );

img_tps1_b_os3_dv2_p1 = reconPwlsTruePathSeekingIndep( sinoSparse, weightsSparse, geomSparse, 'isotv', delta, os, img_sparse_tv_2, img_sparse_tv_1, dv, p, noframes);

figure; imdisp( img_tps1_b_os3_dv2_p1, [0.1 0.3]  );

img_tps2_b_os3_dv2_p1 = reconPwlsTruePathSeekingSelf( sinoSparse, weightsSparse, geomSparse, 'isotv', delta, os, img_sparse_tv_2, img_sparse_tv_1, dv, p, noframes);

figure; imdisp( img_tps2_b_os3_dv2_p1, [0.1 0.3]  );
