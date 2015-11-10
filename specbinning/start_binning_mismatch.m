function start_binning_mismatch
% Start function to run spectrum mismatch simulation. 
% Qiao Yang 
% 2013.10

load 'temp.mat';
close all;
clc;

if exist('results','dir') == 0
    mkdir('results');
end

% Load simulation parameters and datas
[phan, map, roi] = loadMaterialsPhantom(p);

[geom ] = loadProjectionGeometryCT( p );


materialsNames = phan.materialNames;


% Using the 120kV-2mmAl spectrum for estimation of binning parameters
p.Spectra.spectrumKeV = 'spectrum-120kv-2mmAl';
spectrum = loadSpectraCT(p, geom, 2e6);
[binning2.photonCounts, binning2.Phis, binning2.Thetas] = spectrumBinning2Decomposition( spectrum, materialsNames{2}, materialsNames{3} );
[binning3.photonCounts, binning3.Phis, binning3.Thetas] = spectrumBinning3Decomposition( spectrum,materialsNames{2}, materialsNames{3} );


% Start simulations using different spectra
spectrumNames = {'spectrum-80kv-2mmAl';'spectrum-100kv-2mmAl';'spectrum-140kv-2mmAl';...
                 'spectrum-120kv-1mmAl';'spectrum-120kv-2mmAl-1mmTi';'spectrum-120kv-2mmAl-0.5mmCu'};
%celldata = cellstr(data);
for i = 1:1%length(spectrumNames)
    fprintf('\n Start simulating using energy spectrum : "%s.txt" \n\n', spectrumNames{i} );
    p.Spectra.spectrumKeV = spectrumNames{i};
    bhc(map, roi, sinosDirKeV, p, geom, phan,materialsNames, binning2, binning3)
end


%Reconstructions for certain spectrum.
function bhc(map, roi, sinosDirKeV, p, geom, phan, materialsNames, binning2, binning3)

displayResult = true;
saveToRaw = true;

% Compute spectrum
spectrum = loadSpectraCT(p, geom, 2e6);

% Compute sinogram
sinoPhotonCounts = computePhotonCountSinogram( phan, geom, spectrum, sinosDirKeV );

sinoAtt = computeSinogramAttunation( sinoPhotonCounts, spectrum );

% Compute ground truth
[ imgGtAtt, imgGtHu ] = computeGroundTruthCT(phan, spectrum);

[ phiSoft, thetaSoft, muSoft ] = attenuationCoefficientDecompose( materialsNames{2}, spectrum );
[ phiBone, thetaBone, muBone ] = attenuationCoefficientDecompose( materialsNames{3}, spectrum );

%% FBP reconstrunction

imgAttFBP = reconFBP( sinoAtt, geom);

if displayResult 
    imgHuFBP = imgIntensityCorrectionn( imgAttFBP, imgGtHu, map  );
    compare2GroundTruth( imgHuFBP, imgGtHu, roi, map );
    figure('Name', 'FBP KeV Image', 'NumberTitle', 'off');
    imshow( imgHuFBP', roi.windowHu);
    saveas(gcf, ['results/' p.Spectra.spectrumKeV '_FBPw.eps']); 
end
if saveToRaw
write_data_matrix(['results/' p.Spectra.spectrumKeV '_FBPw.raw'],imgHuFBP');
end


% With JS correction
sinoAttJS = sinoAtt + beamHardeningCorrectionJoseph( imgAttFBP, geom, spectrum, materialsNames{2}, materialsNames{3} );

imgAttFBPJS = reconFBP( sinoAttJS, geom);

u0 = imgAttFBPJS;

if displayResult 
    
    imgHuFBPJS = imgIntensityCorrectionn( imgAttFBPJS, imgGtHu, map  );
    compare2GroundTruth( imgHuFBPJS, imgGtHu, roi, map );
    figure('Name', 'FBP JS KeV Image', 'NumberTitle', 'off');
    imshow( imgHuFBPJS', roi.windowHu);
    saveas(gcf, ['results/' p.Spectra.spectrumKeV '_FBPJS.eps']);
end

if saveToRaw
write_data_matrix(['results/' p.Spectra.spectrumKeV '_FBPJS.raw'],imgHuFBPJS');
end

%% PLWS using KeV sinogram
beta = 5e4;
delta= 1e-3;
itnlim = 1;

l = beamHardeningWarterCorrection(sinoAtt, spectrum );

w = computeWeightsPwls( sinoPhotonCounts, spectrum.backgroundEvents );

imgAttPWLS = reconPwlsSeNesterovNusqs( l, w, geom, beta, 'huber', itnlim, delta );

if displayResult
    
    imgHuPWLS =  imgIntensityCorrectionn( imgAttPWLS, imgGtHu, map );
    compare2GroundTruth( imgHuPWLS, imgGtHu, roi, map );
    figure('Name', 'PWLS SE KeV Image', 'NumberTitle', 'off');
    imshow( imgHuPWLS', roi.windowHu);
    saveas(gcf, ['results/' p.Spectra.spectrumKeV '_PWLS.eps']);
    
end

if saveToRaw
write_data_matrix(['results/' p.Spectra.spectrumKeV '_PWLS.raw'],imgHuPWLS');
end

%% PSR using KeV sinogram

%naive soft tissue and bone segementation
[muSoftEff, muBoneEff] = getEffectiveAttuation( u0, muSoft, muBone );
[fsoft, fbone, rho] = materialSegmentationCont( u0 , muSoftEff , muBoneEff );

%Read in segmentation results from K-means for multi-material .
usingKmeansResults = false;
if usingKmeansResults 
segmentation = read_data_matrix('seg-vol.raw', 360, 360);
fMetal = (segmentation == 4 ) ==1;
fsoft = (segmentation == 2 ) == 1;
fbone =(segmentation == 3 ) == 1;
end 


beta = 5e4;
delta= 1e-3;
nitn = 1;

imgAttPSR = reconPolychromaticMaximumLikelihoodSegmLookupTable( sinoPhotonCounts, ...
    spectrum, geom, rho, fsoft, fbone, beta, 'huber', nitn, delta, materialsNames{2}, materialsNames{3});

if displayResult
    imgHuPSR =  imgIntensityCorrectionn( imgAttPSR, imgGtHu, map );
    compare2GroundTruth( imgHuPSR, imgGtHu, roi, map );
    figure('Name', 'PSR full KeV Image', 'NumberTitle', 'off');
    imshow(imgHuPSR', roi.windowHu);
    saveas(gcf, ['results/' p.Spectra.spectrumKeV  '_PSRfull.eps']);
end

if saveToRaw
write_data_matrix(['results/' p.Spectra.spectrumKeV '_PSRfull.raw'],imgHuPSR');
end

    

%% Physice-based Spectrum Binning 2

[muSoftEff, muBoneEff] = getEffectiveAttuation( u0, muSoft, muBone );
[fsoft, fbone, rho] = materialSegmentationCont( u0, muSoftEff , muBoneEff );

fphi = fsoft * phiSoft + fbone * phiBone;
ftheta = fsoft * thetaSoft + fbone * thetaBone;
rho = u0 ./ (  fphi + ftheta );

imgAttD2 = reconPolychromaticMaximumLikelihoodDecompSpectrumBinning( sinoPhotonCounts, ...
    binning2.photonCounts, binning2.Phis, binning2.Thetas, geom, rho, fphi, ftheta, beta, 'huber', nitn, delta );

if displayResult
    imgHuD2 =  imgIntensityCorrectionn( imgAttD2, imgGtHu, map );
    compare2GroundTruth( imgHuD2, imgGtHu, roi, map );
    figure('Name', 'PSR D2 Image', 'NumberTitle', 'off');
    imshow( imgHuD2', roi.windowHu);
    saveas(gcf, ['results/' p.Spectra.spectrumKeV '_PSR_D2_mismatch.eps']);
end

if saveToRaw
write_data_matrix(['results/' p.Spectra.spectrumKeV '_PSR_D2_mismatch.raw'],imgHuD2');
end

%% Physice-based Spectrum Binning 3 

imgAttD3 = reconPolychromaticMaximumLikelihoodDecompSpectrumBinning( sinoPhotonCounts, ...
    binning3.photonCounts, binning3.Phis, binning3.Thetas, geom, rho, fphi, ftheta, beta, 'huber', nitn, delta );

if displayResult
    
    imgHuD3 =  imgIntensityCorrectionn( imgAttD3, imgGtHu, map );
    compare2GroundTruth( imgHuD3, imgGtHu, roi, map );
    figure('Name', 'PSR D3 Image', 'NumberTitle', 'off');
    imshow( imgHuD3', roi.windowHu);
    saveas(gcf, ['results/' p.Spectra.spectrumKeV '_PSR_D3_mismatch.eps']);
    
end

if saveToRaw
write_data_matrix(['results/' p.Spectra.spectrumKeV '_PSR_D3_mismatch.raw'],imgHuD3');
end


%% plot the profile
figure;
plot( [ imgHuPWLS(:,end/2) imgHuPSR(:,end/2) imgHuD2(:,end/2) imgHuD3(:,end/2);]  );
legend('PWLS', 'PSR full', 'PSR 2bins', 'PSR 3bins' );

