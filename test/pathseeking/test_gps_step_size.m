load 'temp.mat';

% Load simulation parameters and datas
[phan, map] = loadXCATPhantom(p);

[ geom ] = loadProjectionGeometryCT( p );

spectrum = loadSpectraCT(p, geom, 1e6);

% compute sinogram

[sinoRaw, tubeCurrentProfile] = simulateCTRawData( phan, geom, spectrum, sinosDir );

sinoAtt = processCTRawData( sinoRaw, spectrum, tubeCurrentProfile);

%Compute ground trut
[imgGtAtt, imgGtHu ] = computeGroundTruth(phan, spectrum);

weights = computeWeightsPwls( sinoRaw, 0, spectrum.electronicNoise );


%% parameters

load 'pwls_image_small_betas.mat'
delta = 1e-3;
numos = 16;


%%

dv = 2e-4;
p = 0.1;

x1_f_1_10 = reconPwlsPathSeeking( sinoAtt, weights, geom, 'huber', delta, numos, img_pwls1, img_pwls2, dv, p);
%figure; imdisp( v1 )

x1_b_1_10 = reconPwlsPathSeeking( sinoAtt, weights, geom, 'huber', delta, numos, img_pwls2, img_pwls1, dv, p );
%figure; imdisp( v2 )

x2_f_1_10 = reconPwlsPathSeeking2( sinoAtt, weights, geom, 'huber', delta, numos, img_pwls1, img_pwls2, dv, p);



%%
dv = 2e-4;
p = 0.2;

x1_f_1_20 = reconPwlsPathSeeking( sinoAtt, weights, geom, 'huber', delta, numos, img_pwls1, img_pwls2, dv, p);
%figure; imdisp( v1 )

x1_b_1_20 = reconPwlsPathSeeking( sinoAtt, weights, geom, 'huber', delta, numos, img_pwls2, img_pwls1, dv, p );
%figure; imdisp( v2 )

x2_f_1_20 = reconPwlsPathSeeking2( sinoAtt, weights, geom, 'huber', delta, numos, img_pwls1, img_pwls2, dv, p);


%%
dv = 2e-4;
p = 0.3;

x1_f_1_30 = reconPwlsPathSeeking( sinoAtt, weights, geom, 'huber', delta, numos, img_pwls1, img_pwls2, dv, p);
%figure; imdisp( v1 )

x1_b_1_30 = reconPwlsPathSeeking( sinoAtt, weights, geom, 'huber', delta, numos, img_pwls2, img_pwls1, dv, p );
%figure; imdisp( v2 )

x2_f_1_30 = reconPwlsPathSeeking2( sinoAtt, weights, geom, 'huber', delta, numos, img_pwls1, img_pwls2, dv, p);

%%
dv = 4e-4;
p = 0.05;

x1_f_2_05 = reconPwlsPathSeeking( sinoAtt, weights, geom, 'huber', delta, numos, img_pwls1, img_pwls2, dv, p);
%figure; imdisp( v1 )

x1_b_2_05 = reconPwlsPathSeeking( sinoAtt, weights, geom, 'huber', delta, numos, img_pwls2, img_pwls1, dv, p );
%figure; imdisp( v2 )

x2_f_2_05 = reconPwlsPathSeeking2( sinoAtt, weights, geom, 'huber', delta, numos, img_pwls1, img_pwls2, dv, p);




%%
dv = 4e-4;
p = 0.1;

x1_f_2_10 = reconPwlsPathSeeking( sinoAtt, weights, geom, 'huber', delta, numos, img_pwls1, img_pwls2, dv, p);
%figure; imdisp( v1 )

x1_b_2_10 = reconPwlsPathSeeking( sinoAtt, weights, geom, 'huber', delta, numos, img_pwls2, img_pwls1, dv, p );
%figure; imdisp( v2 )

x2_f_2_10 = reconPwlsPathSeeking2( sinoAtt, weights, geom, 'huber', delta, numos, img_pwls1, img_pwls2, dv, p);

%%
dv = 4e-4;
p = 0.2;

x1_f_2_20 = reconPwlsPathSeeking( sinoAtt, weights, geom, 'huber', delta, numos, img_pwls1, img_pwls2, dv, p);
%figure; imdisp( v1 )

x1_b_2_20 = reconPwlsPathSeeking( sinoAtt, weights, geom, 'huber', delta, numos, img_pwls2, img_pwls1, dv, p );
%figure; imdisp( v2 )

x2_f_2_20 = reconPwlsPathSeeking2( sinoAtt, weights, geom, 'huber', delta, numos, img_pwls1, img_pwls2, dv, p);



