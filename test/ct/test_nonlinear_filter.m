%% Compare nonlinear filter for FBP denoising
% Meng Wu at Stanford University and University of Erlangen and Nuremberg 
% 2014

%% Load parameters 

% Load current CT system configuration
load 'temp.mat';

% Load numerical phantom parameters
[phan, map] = loadXCATPhantom(p);

% geometry info
geom = loadProjectionGeometryCT( p );

% Geometry info
spectrum = loadSpectraCT(p, geom, 1e6);

% Compute ground truth
[imgGtAtt, imgGtHu ] = computeGroundTruth(phan, spectrum);

%% reconstruction method

recon = @( sino )convertMonoAttToHu( reconFBP( sino, geom, 'ram-lak', 1, 1), spectrum);
reconVar = @( sino, type )convertVarAttToHu( reconFBPVar( sino, geom, 'ram-lak', type, 1), spectrum);
roistd = @( x )std( reshape( x( 221: 260, 201: 240, 2 ), [1600 1] ) );


%% Simulate the raw data

[sinoRaw, tubeCurrentProfile] = simulateCTRawData( phan, geom, spectrum, sinosDir );
sinoRawNoiseless = simulateCTRawData( phan, geom, spectrum, sinosDir, false );

% figure; plot( tubeCurrentProfile );

%% Process the raw data.
% Flat filed normalization and log.

sinoAtt = processCTRawData( sinoRaw, spectrum, tubeCurrentProfile);
sinoAttNoiseless = processCTRawData( sinoRawNoiseless, spectrum, tubeCurrentProfile);

%% FBP reconstruction

imgFBP = recon( sinoAtt );
figure; imshow( imgFBP(:, :, 2), [-150 150] ); 
title( [ 'FBP reconstruction ' num2str( roistd( imgFBP ) ) ' HU' ]);

imgFBPNoiseless = recon( sinoAttNoiseless );
figure; imshow( imgFBPNoiseless(:, :, 2), [-150 150] ); 
title( [ 'FBP noiseless reconstruction ' num2str( roistd( imgFBPNoiseless ) ) ' HU' ]);

figure; imshow( imgFBP(:, :, 2) - imgFBPNoiseless(:, :, 2), [-50 50] );
title( 'Noise in FBP' );

%% Reconstruct variance map
% L. Zhu and J. StarLack, “A Practical Reconstruction Algorithm for CT
% Noise Variance Maps Using FBP Reconstruction,” in SPIE Medical Imaging,
% 2007, vol. 6510, pp. 651023–651023–8.    

sinoVar = computeProjectionVariance( sinoRaw, spectrum.electronicNoise, spectrum.detectorGain );
% figure; imdisp( sinoVar(3 ,:,:), [0 1e-3 ]);

imgVarH = reconVar( sinoVar, 'horiz') ;
figure; imshow( imgVarH(:,:,2), [0 1000] );
title 'Variance from horizontal rays';

imgVarV = reconVar( sinoVar, 'vert');
figure; imshow( imgVarV(:,:,2), [0 1000] );
title 'Variance from vertical rays';

imgVarAll = imgVarH + imgVarV;
figure; imshow( imgVarAll(:,:,2), [0 1000] );
title 'Variance from all rays';

% roiVar = squeeze( imgVarAll( 200: 240, 220: 260, 2 ) );
% mean(  sqrt( roiVar(:) ) )

%% Projction Space Bilateral Filter
% A. Manduca, L. Yu, J. D. Trzasko, N. Khaylova, J. M. Kofler, C. M.
% McCollough, and J. G. Fletcher, “Projection space denoising with
% bilateral filtering and CT noise modeling for dose reduction in CT,” Med.
% Phys., vol. 36, no. 11, p. 4911, Nov. 2009.    

sinoAttBilt = bilateralFilterProjectionDomain( sinoAtt, geom, 2, 1, 1 );
imgProjBilt = recon( sinoAttBilt );

figure; imshow( imgProjBilt(:, :, 2), [-150 150] );
title( [ 'Projection Domian Bileteral Filter ' num2str( roistd( imgProjBilt ) ) ' HU' ]);

figure; imshow( imgProjBilt(:, :, 2) - imgFBPNoiseless(:, :, 2), [-50 50] );
title( 'Noise in Projection Domian Bileteral Filter' );


%% Bilateral Filter in Image Domain 
% J. C. R. Giraldo, Z. S. Kelm, L. S. Guimaraes, L. Yu, J. G. Fletcher,
% B. J. Erickson, and C. H. McCollough, “Comparative study of two image
% space noise reduction methods for computed tomography: bilateral filter
% and nonlocal means.,” IEEE Eng. Med. Biol. Soc. Conf., vol. 2009, pp.
% 3529–32, Jan. 2009.    
%
% A. R. Al-Hinnawi, M. Daear, and S. Huwaijah, “Assessment of bilateral
% filter on 1/2-dose chest-pelvis CT views.,” Radiol. Phys. Technol., vol.
% 6, no. 2, pp. 385–98, Jul. 2013.    

imgBilt = bilateralFilterImageDomain( imgFBP, 2, 1, 50 );

figure; imshow( imgBilt(:, :, 2), [-150 150] ); 
title( [ 'Image Domian Bileteral Filter ' num2str( roistd( imgBilt ) ) ' HU' ]);

figure; imshow( imgBilt(:, :, 2) - imgFBPNoiseless(:, :, 2), [-50 50] );
title( 'Noise in Image Domian Bileteral Filter' );

%% Anisotrpic adaptive bilateral filter in image domian
%
% By Meng Wu

imgAnisoBilt = anisotropicAdaptiveBilateralFilterImageDomain( imgFBP, imgVarH, imgVarV, 2, 1, 75, 1 );

figure; imshow( imgAnisoBilt(:, :, 2), [-150 150] ); 
title( [ 'Image Domian Anisotropic Bileteral Filter ' num2str( roistd( imgAnisoBilt ) ) ' HU' ]);

figure; imshow( imgAnisoBilt(:, :, 2) - imgFBPNoiseless(:, :, 2), [-50 50] );
title( 'Noise in Image Domian Anisotropic Bileteral Filter' );

%% Nonlocal Mean Filter in Image Domain
% J. C. R. Giraldo, Z. S. Kelm, L. S. Guimaraes, L. Yu, J. G. Fletcher,
% B. J. Erickson, and C. H. McCollough, “Comparative study of two image
% space noise reduction methods for computed tomography: bilateral filter
% and nonlocal means.,” IEEE Eng. Med. Biol. Soc. Conf., vol. 2009, pp.
% 3529–32, Jan. 2009.    

imgNLM = nonLocalMeanFilterImageDomain( imgFBP, 2, 1, 1, 20 );

figure; imshow( imgNLM(:, :, 2), [-150 150] ); 
title( [ 'NonLocal Mean Filter ' num2str( roistd( imgNLM ) ) ' HU' ]);

figure; imshow( imgNLM(:, :, 2) - imgFBPNoiseless(:, :, 2), [-50 50] );
title( 'Noise in NonLocal Mean Filter' );

%% Struture Adaptive Sinogram Filter 
% Based on M. Balda, J. Hornegger, and B. Heismann, “Ray contribution masks for
% structure adaptive sinogram filtering.,” IEEE Trans. Med. Imaging, vol.
% 31, no. 6, pp. 1228–39, Jun. 2012.  

% sinoAttSAS =  structurAdaptiveSinogramFilter( sinoAtt, geom, 3, 1 );

imgFBPSAS = recon( sinoAttSAS );

figure; imshow( imgFBPSAS(:, :, 2), [-150 150] );
title( [ 'Struture adaptive sinogram filter ' num2str( roistd( imgFBPSAS ) ) ' HU' ]);

figure; imshow( imgFBPSAS(:, :, 2) - imgFBPNoiseless(:, :, 2), [-50 50] );
title( 'Noise in NonLocal Mean Filter' );

return;
%% Local Quadratic Programming Filter
%
% By Meng Wu

% imgFBPLQP = localQuadFilterImageDomain( imgFBP, imgVarAll, 3, 1 );

figure; imshow( imgFBPLQP(:, :, 2), [-150 150] ); 
title( [ 'Local Quadratic Programming Filter ' num2str( roistd( imgFBPLQP ) ) ' HU' ]);


figure; imshow( imgFBPLQP(:, :, 2) - imgFBPNoiseless(:, :, 2), [-50 50] );
title( 'Noise in Local Quadratic Programming Filter' );



