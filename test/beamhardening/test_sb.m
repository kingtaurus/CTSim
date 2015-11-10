load 'temp.mat';
displayResult = true;

% Load simulation parameters and datas

[phan, map, roi] = loadMaterialsPhantom(p);

[geom ] = loadProjectionGeometryCT( p );

spectrum = loadSpectraCT(p, geom, 3e6);

% compute sinogram
sinoPhotonCounts = simulatePhotonCountingData( phan, geom, spectrum, sinosDirKeV );
sinoAtt = computePhotonCountingSinogramAttunation( sinoPhotonCounts, spectrum );

% Compute ground trut
[ imgGtAtt, imgGtHu ] = computeGroundTruthCT(phan, spectrum);

sinoAttWC = beamHardeningWarterCorrection(sinoAtt, spectrum );


material1 = phan.materialNames{2};
material2 = phan.materialNames{3};
[ phiSoft, thetaSoft, muSoft, Phi, Theta ] = attenuationCoefficientDecompose( material1, spectrum );
[ phiBone, thetaBone, muBone ] = attenuationCoefficientDecompose( material2, spectrum );

%% FBP reconstrunction with water correction

imgAttFBP = reconFBP( sinoAttWC, geom);
imgHuFBP_WC2 = imgIntensityCorrectionn( imgAttFBP, imgGtHu, map  );
compare2GroundTruth( imgHuFBP_WC2, imgGtHu, roi, map, 'fbp_wc' );


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




%% PSR

beta = 5e4;
delta= 1e-3;
nitn = 50;

[muSoftEff, muBoneEff] = getEffectiveAttuation( u0, muSoft, muBone );
[fsoft, fbone] = materialSegmentationLinear( u0 , muSoftEff , muBoneEff );
fphi = fsoft * phiSoft + fbone * phiBone;
ftheta = fsoft * thetaSoft + fbone * thetaBone;
rho = u0 ./ (  fphi + ftheta );

imgAttPSR = reconPolychromaticMaximumLikelihoodSegmLookupTable( sinoPhotonCounts, ...
    spectrum, geom, rho, fsoft, fbone, beta, 'huber', nitn, delta, material1, material2);

imgHuPSR =  imgIntensityCorrectionn( imgAttPSR, imgGtHu, map );
compare2GroundTruth( imgHuPSR, imgGtHu, roi, map, 'psr_full' );



%% Generate a range of material thicknesses

[phis, thetas, thickness1, thickness2] = generateMaterialThicknessMap( 32, spectrum, material1, material2,  40, 10 );
projThicknessCounts = computeEffectiveThickness( fsoft, fbone, rho, geom, spectrum, material1, material2, phis, thetas );
% weights for spectrum binning
weights = projThicknessCounts / sum(projThicknessCounts) * length( projThicknessCounts );


%% Spectrum binning with 2 energy bins

[Bins_tsb2, Phis_tsb2, Thetas_tsb2] = thresholdBasedSpectrumBinning2( spectrum, 1, phis, thetas );
[Bins_gsb2, Phis_gsb2, Thetas_gsb2] = generalizedSpectrumBinning( spectrum,  1, 2, phis, thetas );
[Bins_wgsb2, Phis_wgsb2, Thetas_wgsb2] = generalizedSpectrumBinning( spectrum,  weights, 2, phis, thetas );


%% Spectrum binning with 3 energy bins

[Bins_tsb3, Phis_tsb3, Thetas_tsb3] = thresholdBasedSpectrumBinning3( spectrum, 1, phis, thetas );
[Bins_wtsb3, Phis_wtsb3, Thetas_wtsb3] = thresholdBasedSpectrumBinning3( spectrum, weights, phis, thetas );
[Bins_gsb3, Phis_gsb3, Thetas_gsb3] = generalizedSpectrumBinning( spectrum,  1, 3, phis, thetas );
[Bins_wgsb3, Phis_wgsb3, Thetas_wgsb3] = generalizedSpectrumBinning( spectrum,  weights, 3, phis, thetas );

%%
[Bins_wgsb4, Phis_wgsb4, Thetas_wgsb4] = generalizedSpectrumBinning( spectrum,  weights, 4, phis, thetas );
showBinnigErrorMap( spectrum, Bins_wgsb4, Phis_wgsb4, Thetas_wgsb4, [-0.1 0.3]  );

%% bowtie info

if spectrum.useBowtie
    bowtieThickness = spectrum.bowtieThickness / 10;
    [ phiF, thetaF ] = attenuationCoefficientDecompose( spectrum.bowtieMaterial, spectrum );
    
    phisFiltraion = bowtieThickness * phiF;
    thetaFiltrations = bowtieThickness * thetaF;
    
else
    phisFiltraion = 0;
    thetaFiltrations = 0;
    
end

%%
beta = 1e5;
delta= 1e-3;
nitn = 50;
numos = 8;

imgAttSB = reconPsrSpectrumBinningNesterovSqs( sinoPhotonCounts, ...
    Bins_wgsb4, Phis_wgsb4, Thetas_wgsb4, geom, rho, fphi, ftheta, phisFiltraion, thetaFiltrations, ...
    beta, 'huber', nitn, delta, numos );

imgHuSB =  imgIntensityCorrectionn( imgAttSB, imgGtHu, map );
compare2GroundTruth( imgHuSB, imgGtHu, roi, map, 'pml_sb' );

%%
beta = 1e5;
delta= 1e-3;
nitn = 50;
numos = 8;

imgAttSB = reconPsrSpectrumBinningADMMSqs( sinoPhotonCounts, ...
    Bins_wgsb4, Phis_wgsb4, Thetas_wgsb4, geom, rho, fphi, ftheta, phisFiltraion, thetaFiltrations, ...
    beta, 'huber', nitn, delta, numos );

imgHuSB =  imgIntensityCorrectionn( imgAttSB, imgGtHu, map );
compare2GroundTruth( imgHuSB, imgGtHu, roi, map, 'pml_sb' );


%% plot the profile

figure;
plot( [ imgHuFBP_JS(:,end/2,end/2) imgHuPWLS(:,end/2,end/2) imgHuPSR(:,end/2,end/2) imgHuSB(:,end/2,end/2) ], 'linewidth', 1 );
legend( 'FBP JS',  'PWLS', 'PSR-full', 'PML-SB',  'Location', 'SouthEast' );
ylabel( 'HU', 'fontSize',14);
saveas(gcf, 'Profiles.eps','epsc');



