clc;
load 'temp.mat';

dataFileName = 'MDCT_Liver_120kVp_200mAs_FS3_DCS5_VS3_T3s_m';

noiselessFileName = 'MDCT_Liver_120kVp_200mAs_FS3_DCS5_VS3_Noiseless';

[geom ] = loadProjectionGeometryCT( p );

[sinoAtt, sinoPC]= loadLiverCTSinogram( dataFileName, geom );

weights = computeWeightsPwls( sinoPC, 0, 0);

%% FBP reconstruction

imgAttFBP = reconFBP( sinoAtt, geom, 'hamming');

imgFBP = (( imgAttFBP(:,:,2) / 0.199 * 1000 ) - 1000 );
figure; imshow( imgFBP', [-150 250] ); colormap gray;
title 'FBP';


%%
sinoVar = computeProjectionVariance( sinoPC / 10 , 0 );
imgVar = reconEquiangularFDKvarianceMap( sinoVar, geom, 'hamming');

imgVar = imgVar(:,:,2);
imgVar = imgVar * ( 1 / 0.199 * 1000 )^2;


imdisp( imgVar); colormap gray;
title 'Variance map';


%% FBP reconstruction

sinoAttBL = bilateralFilterProjectionDomain( sinoAtt, sinoPC, geom, 5, 2, 0.1 );

imgAttBL = reconFBP( sinoAttBL, geom, 'hamming');

imgBL = (( imgAttBL(:,:,2) / 0.199 * 1000 ) - 1000 );
figure; imshow( imgBL', [-150 250] ); colormap gray;
title 'FBP BLP';

%% FBP reconstruction

imgBL = bilateralFilterImageDomain( imgFBP, 5, 2, 100 );

imdisp( imgBL, [-150 250] ); 

%% FBP reconstruction

imgAttNLM = NLMeansFilterImageDomain( imgAttFBP, 7, 2, 0.02, 2 );

imgNLM = (( imgAttNLM(:,:,2) / 0.199 * 1000 ) - 1000 );
figure; imshow( imgNLM', [-150 250] ); colormap gray;
title 'FBP NLM';


%% PLWS using KeV sinogram
nitn = 30;
beta = 1e5;
delta = 1e-3;

imgAttPWLS = reconPwlsSeNesterovNusqs( sinoAtt, weights, geom, beta, 'huber', nitn, delta, false);


imgPWLS = (( imgAttPWLS(:,:,2) / 0.199 * 1000 ) - 1000 );
figure; imshow( imgPWLS', [-150 250] ); colormap gray;
title 'PWLS';


%% PLWS using KeV sinogram
nitn = 30;
beta = 1e4;
delta = 1e-3;

imgAttPWLS = reconPwlsSeNesterovNusqs( sinoAtt, weights, geom, beta, 'huber', nitn, delta, false);


imgPWLS = (( imgAttPWLS(:,:,2) / 0.199 * 1000 ) - 1000 );
figure; imshow( imgPWLS', [-150 250] ); colormap gray;
title 'PWLS';

return;

