%% test_pathseeking_nps
% Compare several reconstruction methods for sparse view CT
% Meng Wu at University of Erlangen-Nuremburg
% 2014.10

%% Simulate normal CT data
clear;
close all;
load 'temp.mat';


% Load simulation parameters and datas
[phan, map] = loadMaterialsDensityPhantom(p);

% geom = loadProjectionGeometryShortScan( p, 223, 90 );

geom = loadProjectionGeometryCT( p );

spectrum = loadSpectraCT(p, geom, 2e5);

for id = 1 
    
    % compute sinogram
    [sinoRaw, tubeCurrentProfile] = simulateCTRawData( phan, geom, spectrum, sinosDir, 1, id, 2 );
    
    sinoAtt = processCTRawData( sinoRaw, spectrum, tubeCurrentProfile);
    
    weights = computeWeightsPwls( sinoRaw, 0, spectrum.electronicNoise );
    
    img_fbp_sharp = reconFBP( sinoAtt, geom, 'ram-lak', 1, 1);
    %figure; imshow( img_fbp_sharp(:,:,ceil(end/2)), map.windowAtt );
    %title 'FBP reconstruction with ram-lak kernel';
    
    img_fbp_soft = reconFBP( sinoAtt, geom, 'hamming', 1, 1);
    %figure; imshow( img_fbp_soft(:,:,ceil(end/2)), map.windowAtt );
    %title 'FBP reconstruction with hamming window';
    
    img_fbp_soft1 = reconFBP( sinoAtt, geom, 'hamming', 1,  0.5 );
    %figure; imshow( img_fbp_soft(:,:,ceil(end/2)), map.windowAtt );
    %title 'FBP reconstruction with hamming window';
    
    %%

    nitn = 20;
    numos =16;
    beta1 = 5e3;
    beta2 = 1e5;
    beta0 = 2e4;
    delta = 1e-3;
    
    img_pwls_beta1 = reconPwlsLALMOs14( sinoAtt, weights, geom, beta1, 'huber', nitn, delta, numos, img_fbp_soft  );
    figure; imdisp( img_pwls_beta1, map.windowAtt );
    %title( ['PWLS reconstruction beta = ' num2str(beta1)] );
    
    img_pwls_beta2 = reconPwlsLALMOs14( sinoAtt, weights, geom, beta2, 'huber', nitn, delta, numos, img_fbp_soft  );
    figure; imdisp( img_pwls_beta2, map.windowAtt );
    %title( ['PWLS reconstruction beta = ' num2str(beta2)] );
    
    img_pwls_beta0 = reconPwlsLALMOs14( sinoAtt, weights, geom, beta0, 'huber', nitn, delta, numos, img_fbp_soft  );
    figure; imdisp( img_pwls_beta0, map.windowAtt );
    %title( ['PWLS reconstruction beta = ' num2str(beta0)] );
    
    %%
        return;
    
        
    dv = 4e-4;
    p = 0.2;
    os = 4;
    
    img_aps = reconPwlsApproxPathSeeking( sinoAtt, weights, geom, 'huber', delta, os, img_pwls_beta1, img_pwls_beta2, dv, p, 40);
    
    img_tps = reconPwlsTruePathSeekingSelf( sinoAtt, weights, geom, 'huber', delta, os, img_pwls_beta1, img_pwls_beta2, dv, p, 40);
    
    save( [sinosDir 'water_' num2str(id) '.mat'], 'img_fbp_sharp', 'img_fbp_soft', 'img_pwls_beta1', 'img_pwls_beta2' , 'img_pwls_beta0', 'img_aps', 'img_tps' );
    
end