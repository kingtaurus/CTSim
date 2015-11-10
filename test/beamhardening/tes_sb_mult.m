load 'temp.mat';
displayResult = true;
close all;
clc;

% Load simulation parameters and datas
[phan, map, roi] = loadMaterialsPhantom(p);

[geom ] = loadProjectionGeometryCT( p );

spectrum = loadSpectraCT(p, geom, 1e6);

% compute sinogram
sinoPhotonCounts = simulatePhotonCountingData( phan, geom, spectrum, sinosDirKeV );
sinoAtt = computePhotonCountingSinogramAttunation( sinoPhotonCounts, spectrum );

% Compute ground trut
[ imgGtAtt, imgGtHu ] = computeGroundTruthCT(phan, spectrum);
sinoAttWC = beamHardeningWarterCorrection(sinoAtt, spectrum );

[ phiSoft, thetaSoft, muSoft, Phi, Theta ] = attenuationCoefficientDecompose( 'Tissue_Soft_ICRU-44', spectrum );
[ phiBone1, thetaBone1, muBone1 ] = attenuationCoefficientDecompose( 'B-100_Bone-Equivalent_Plastic', spectrum );
[ phiBone2, thetaBone2, muBone2 ] = attenuationCoefficientDecompose( 'Bone_Cortical_ICRU-44', spectrum );
[ phiAl, thetaAl, muAl ] = attenuationCoefficientDecompose( 'Aluminum', spectrum );
[ phiCa, thetaCa, muCa ] = attenuationCoefficientDecompose( 'Ca', spectrum );



material1 = phan.materialNames{2};
material2 = phan.materialNames{3};

% projThicknessCounts = computeEffectiveThickness( sinoAtt, geom, spectrum, 16, material1, material2  );
% weights = projThicknessCounts / sum(projThicknessCounts) * length( projThicknessCounts );

%% FBP reconstrunction


imgAttFBP = reconFBP( sinoAttWC, geom);
imgHuFBP_WC = imgIntensityCorrectionn( imgAttFBP, imgGtHu, map  );
compare2GroundTruth( imgHuFBP_WC, imgGtHu, roi, map, 'fbp_wc' );


% FBP reconstrunction with JS correction
sinoAttJS = sinoAtt + beamHardeningCorrectionJoseph( imgAttFBP, geom, spectrum, material1, material2 );
imgAttFBP = reconFBP( sinoAttJS, geom);
imgHuFBP_JS = imgIntensityCorrectionn( imgAttFBP, imgGtHu, map  );
compare2GroundTruth( imgHuFBP_JS, imgGtHu, roi, map, 'fbp_js' );


u0 = imgAttFBP;


%% pwls
beta = 5e4;
delta= 1e-3;
nitn = 30;
l = beamHardeningWarterCorrection(sinoAtt, spectrum );
w = computeWeightsPwls( sinoPhotonCounts, 0, spectrum.backgroundEvents );

imgAttPWLS = reconPwlsSeNesterovSqs( l, w, geom, beta, 'huber', nitn, delta, 8 );

imgHuPWLS =  imgIntensityCorrectionn( imgAttPWLS, imgGtHu, map );
compare2GroundTruth( imgHuPWLS, imgGtHu, roi, map, 'pwls' );

%%
beta = 5e4;
delta= 1e-3;
nitn = 50;

[ phiSoft, thetaSoft, muSoft, Phi, Theta ] = attenuationCoefficientDecompose( material1, spectrum );
[ phiBone, thetaBone, muBone ] = attenuationCoefficientDecompose( material2, spectrum );
[muSoftEff, muBoneEff] = getEffectiveAttuation( u0, muSoft, muBone );
[fsoft, fbone] = materialSegmentationLinear( u0 , muSoftEff , muBoneEff );
fphi = fsoft * phiSoft + fbone * phiBone;
ftheta = fsoft * thetaSoft + fbone * thetaBone;
rho = u0 ./ (  fphi + ftheta );

imgAttPSR = reconPolychromaticMaximumLikelihoodSegmLookupTable( sinoPhotonCounts, ...
    spectrum, geom, rho, fsoft, fbone, beta, 'huber', nitn, delta, material1, material2);

imgHuPSR =  imgIntensityCorrectionn( imgAttPSR, imgGtHu, map );
compare2GroundTruth( imgHuPSR, imgGtHu, roi, map, 'psr_full' );


%% PSR
beta = 2e4;
delta= 1e-3;
nitn = 50;

%[muSoftEff, muBoneEff] = getEffectiveAttuation( u0, muSoft, muBone );
[fsoft, fbone, rho] = materialSegmentationCont( u0 , 0.23 , 0.56 );

imgAtt = reconPolychromaticMaximumLikelihoodSegmLookupTable( sinoPhotonCounts, ...
    spectrum, geom, rho, fsoft, fbone, beta, 'huber', nitn, delta, material1, material2);

if displayResult
    imgHu =  imgIntensityCorrectionn( imgAtt, imgGtHu, map );
    compare2GroundTruth( imgHu, imgGtHu, roi, map );
    figure('Name', 'PSR full', 'NumberTitle', 'off');
    imshow( imgHu', roi.windowHu);
    saveas(gcf, 'PSR_full.eps');
    img_PSR = imgHu;
end

%%

%[muSoftEff, muBoneEff] = getEffectiveAttuation( u0, muSoft, muBone );
%[fsoft, fbone, rho] = materialSegmentationCont( u0 , muSoftEff , muBoneEff );
[ fsoft, fbone1, fbone2, fAl, fCa ] =  materialSegmentationMultBHC( u0, muBone1, muBone2, muAl, muCa );

fphi = fsoft * phiSoft + fbone1 * phiBone1 + fbone2  * phiBone2 + fAl * phiAl + fCa * phiCa;
ftheta = fsoft * thetaSoft + fbone1 * thetaBone1 + fbone2  * thetaBone2 + fAl * thetaAl + fCa * thetaCa;
rho = u0 ./ (  fphi + ftheta );
figure; imagesc( fphi );


material1 = 'Tissue_Soft_ICRU-44';
%material2 = 'Bone_Cortical_ICRU-44';
material2 =  'Ca';

%% Spectrum Binning 2
beta = 2e4;
delta= 1e-3;
nitn = 50;


[photonCounts2, Phis2, Thetas2] = spectrumBinning2Decomposition( spectrum, 1,  material1, material2, 40, 10 );
saveas(gcf, 'SB2.eps','epsc');
showBinnigErrorMap( spectrum, photonCounts2, Phis2, Thetas2, [-0.1 0.2]  );
saveas(gcf, 'ErrorD2.eps','epsc');

imgAtt = reconPolychromaticMaximumLikelihoodDecompSpectrumBinning( sinoPhotonCounts, ...
    photonCounts2, Phis2, Thetas2, geom, rho, fphi, ftheta, beta, 'huber', nitn, delta );

if displayResult
    imgHu =  imgIntensityCorrectionn( imgAtt, imgGtHu, map );
    compare2GroundTruth( imgHu, imgGtHu, roi, map );
    figure('Name', 'PSR D2 Image', 'NumberTitle', 'off');
    imshow( imgHu', roi.windowHu);
    saveas(gcf, 'PSRD2.eps');
    img_PSR_D2 = imgHu;
end


% % weighted
% [photonCounts2, Phis2, Thetas2] = spectrumBinning2Decomposition( spectrum, weights,  material1, material2 );
% showBinnigErrorMap( spectrum, photonCounts2, Phis2, Thetas2 );
% saveas(gcf, 'ErrorD2W.eps','epsc');
% 
% imgAtt = reconPolychromaticMaximumLikelihoodDecompSpectrumBinning( sinoPhotonCounts, ...
%     photonCounts2, Phis2, Thetas2, geom, rho, fphi, ftheta, beta, 'huber', nitn, delta );
% 
% if displayResult
%     imgHu =  imgIntensityCorrectionn( imgAtt, imgGtHu, map );
%     compare2GroundTruth( imgHu, imgGtHu, roi, map );
%     figure('Name', 'PSR D2W Image', 'NumberTitle', 'off');
%     imshow( imgHu', roi.windowHu);
%     saveas(gcf, 'PSRD2W.eps');
%     img_PSR_D2W = imgHu;
% end

[photonCounts2, Phis2, Thetas2] = generalizedSpectrumBinning2( spectrum, 1, material1, material2, 40, 10);
showBinnigErrorMap( spectrum, photonCounts2, Phis2, Thetas2, [-0.1 0.2] );
saveas(gcf, 'ErrorG2.eps' ,'epsc');

imgAtt = reconPolychromaticMaximumLikelihoodDecompSpectrumBinning( sinoPhotonCounts, ...
    photonCounts2, Phis2, Thetas2, geom, rho, fphi, ftheta, beta, 'huber', nitn, delta );

if displayResult
    imgHu =  imgIntensityCorrectionn( imgAtt, imgGtHu, map );
    compare2GroundTruth( imgHu, imgGtHu, roi, map );
    figure('Name', 'PSR G2 Image', 'NumberTitle', 'off');
    imshow( imgHu', roi.windowHu);
    saveas(gcf, 'PSRG2.eps');
    img_PSR_G2 = imgHu;
end

% % weighted
% [photonCounts2, Phis2, Thetas2] = generalizedSpectrumBinning2( spectrum, weights, material1, material2);
% showBinnigErrorMap( spectrum, photonCounts2, Phis2, Thetas2 );
% saveas(gcf, 'ErrorG2W.eps','epsc');
% 
% imgAtt = reconPolychromaticMaximumLikelihoodDecompSpectrumBinning( sinoPhotonCounts, ...
%     photonCounts2, Phis2, Thetas2, geom, rho, fphi, ftheta, beta, 'huber', nitn, delta );
% 
% if displayResult
%     imgHu =  imgIntensityCorrectionn( imgAtt, imgGtHu, map );
%     compare2GroundTruth( imgHu, imgGtHu, roi, map );
%     figure('Name', 'PSR G2W Image', 'NumberTitle', 'off');
%     imshow( imgHu', roi.windowHu);
%     saveas(gcf, 'PSRG2W.eps');
%     img_PSR_G2W = imgHu;
% end



%% Spectrum Binning 3
[photonCounts3, Phis3, Thetas3] = spectrumBinning3Decomposition( spectrum, 1, material1, material2, 40, 10);
saveas(gcf, 'SB3.eps','epsc');
showBinnigErrorMap( spectrum, photonCounts3, Phis3, Thetas3, [-0.1 0.2] );
saveas(gcf, 'ErrorD3.eps','epsc');

imgAtt = reconPolychromaticMaximumLikelihoodDecompSpectrumBinning( sinoPhotonCounts, ...
    photonCounts3, Phis3, Thetas3, geom, rho, fphi, ftheta, beta, 'huber', nitn, delta );

if displayResult
    imgHu =  imgIntensityCorrectionn( imgAtt, imgGtHu, map );
    compare2GroundTruth( imgHu, imgGtHu, roi, map );
    figure('Name', 'PSR D3 Image', 'NumberTitle', 'off');
    imshow( imgHu', roi.windowHu);
    saveas(gcf, 'PSRD3.eps');
    img_PSR_D3 = imgHu;
end


% [photonCounts3, Phis3, Thetas3] = spectrumBinning3Decomposition( spectrum, weights, material1, material2 );
% showBinnigErrorMap( spectrum, photonCounts3, Phis3, Thetas3 );
% saveas(gcf, 'ErrorD3W.eps','epsc');
% 
% imgAtt = reconPolychromaticMaximumLikelihoodDecompSpectrumBinning( sinoPhotonCounts, ...
%     photonCounts3, Phis3, Thetas3, geom, rho, fphi, ftheta, beta, 'huber', nitn, delta );
% 
% if displayResult
%     imgHu =  imgIntensityCorrectionn( imgAtt, imgGtHu, map );
%     compare2GroundTruth( imgHu, imgGtHu, roi, map );
%     figure('Name', 'PSR D3W Image', 'NumberTitle', 'off');
%     imshow( imgHu', roi.windowHu);
%     saveas(gcf, 'PSR_D3W.eps');
%     img_PSR_D3W = imgHu;
% end


[photonCounts3, Phis3, Thetas3] = generalizedSpectrumBinning3( spectrum, 1, material1, material2, 40, 10 );
showBinnigErrorMap( spectrum, photonCounts3, Phis3, Thetas3, [-0.1 0.2] );
saveas(gcf, 'ErrorG3.eps','epsc');

imgAtt = reconPolychromaticMaximumLikelihoodDecompSpectrumBinning( sinoPhotonCounts, ...
    photonCounts3, Phis3, Thetas3, geom, rho, fphi, ftheta, beta, 'huber', nitn, delta );

if displayResult
    imgHu =  imgIntensityCorrectionn( imgAtt, imgGtHu, map );
    compare2GroundTruth( imgHu, imgGtHu, roi, map );
    figure('Name', 'PSR G3 Image', 'NumberTitle', 'off');
    imshow( imgHu', roi.windowHu);
    saveas(gcf, 'PSRG3.eps');
    img_PSR_G3 = imgHu;
end

% [photonCounts3, Phis3, Thetas3] = generalizedSpectrumBinning3( spectrum, weights, material1, material2 );
% showBinnigErrorMap( spectrum, photonCounts3, Phis3, Thetas3 );
% saveas(gcf, 'ErrorG3W.eps','epsc');
% 
% imgAtt = reconPolychromaticMaximumLikelihoodDecompSpectrumBinning( sinoPhotonCounts, ...
%     photonCounts3, Phis3, Thetas3, geom, rho, fphi, ftheta, beta, 'huber', nitn, delta );
% 
% if displayResult
%     imgHu =  imgIntensityCorrectionn( imgAtt, imgGtHu, map );
%     compare2GroundTruth( imgHu, imgGtHu, roi, map );
%     figure('Name', 'PSR G3W Image', 'NumberTitle', 'off');
%     imshow( imgHu', roi.windowHu);
%     saveas(gcf, 'PSRG3W.eps');
%     img_PSR_G3W = imgHu;
% end





%% plot the profile

figure;
plot( [ img_FBP_JS(:,end/2) img_PWLS(:,end/2) img_PSR(:,end/2) img_PSR_D2(:,end/2) img_PSR_G2(:,end/2) img_PSR_D3(:,end/2) img_PSR_G3(:,end/2)]  );
legend( 'FBP JS',  'PWLS', 'PSR-full', 'PSR-SB2', 'PSR-GSB2', 'PSR-SB3', 'PSR-GSB3',  'Location', 'SouthEast' );
ylabel( 'HU', 'fontSize',14);
saveas(gcf, 'Profiles.eps','epsc');

% 
