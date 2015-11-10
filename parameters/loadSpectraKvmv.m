function [spectrumKeV, spectrumMeV] = loadSpectraKvmv(p, geomKeV, geomMeV, photonsTotalKeV, photonsTotalMeV )
% Load Spectra
%
% Copyright (c) 2010-2012 by Andreas Keil, Stanford University.
% Modified by Meng Wu at 2012.9

printOut = true;

if nargin < 4
    photonsTotalKeV = 1e6;
    photonsTotalMeV = 5e5;
end

steradiansPerMm2At1m = solidAnglePyramidFromDistanceAndBaseLengths(1, 0.001, 0.001);


%% kV spectrum
[energyBinLabelsKeV, ~, photonsPerEnergyBinKeV, ~, ~, energyAverageKeV] = ...
    readSpectrum([p.Paths.spectraDir p.Spectra.spectrumKeV '.txt'], 'per bin', 1);

photonsPerEnergyBinKeV      = photonsPerEnergyBinKeV * photonsTotalKeV / sum(photonsPerEnergyBinKeV(:));
photonsPerMm2At1mKeV        = photonsTotalKeV * (geomKeV.SDD/1000)^2 / ( geomKeV.detSpacing(1) * geomKeV.detSpacing(end));
photonsPerSteradianKeV      = photonsPerMm2At1mKeV / steradiansPerMm2At1m;

spectrumKeV.energyBinLabels         = energyBinLabelsKeV;
spectrumKeV.photonsPerEnergyBin     = photonsPerEnergyBinKeV ;
spectrumKeV.photonsTotal            = photonsTotalKeV ;
spectrumKeV.energyAverage           = energyAverageKeV;
spectrumKeV.photonsPerMm2At1m       = photonsPerMm2At1mKeV;
spectrumKeV.photonsPerSteradian     = photonsPerSteradianKeV ;
spectrumKeV.DQE                     = p.Detector.detectorConversionEfficiencyKeV;
spectrumKeV.backgroundEvents        = 50;



%% add bowtie filter into kV spectrum

if ~strcmpi(p.Bowtie.shapeType, 'none')
    
    spectrumKeV.useBowtie       = true;
    spectrumKeV.bowtieThickness = loadBowtieFilter( p, geomKeV );
    spectrumKeV.bowtieMaterial  =  p.Bowtie.material;
    spectrumKeV.flatFieldRatio  = zeros( size(spectrumKeV.bowtieThickness) );
    
    if length( size(spectrumKeV.bowtieThickness)) == 1
        for i = 1:length( spectrumKeV.flatFieldRatio)
            spectrumKeV.flatFieldRatio(i) = sum( computeResultingPhotons( photonsPerEnergyBinKeV, ...
                energyBinLabelsKeV, p.Bowtie.material, spectrumKeV.bowtieThickness(i) )) / photonsTotalKeV;
        end
    else
        for i = 1:size( spectrumKeV.flatFieldRatio, 2)
            spectrumKeV.flatFieldRatio(:, i) = sum( computeResultingPhotons( photonsPerEnergyBinKeV, ...
                energyBinLabelsKeV, p.Bowtie.material, spectrumKeV.bowtieThickness(ceil(end/2), i) )) / photonsTotalKeV;
        end
        
    end
    
    spectrumKeV.photonsPerEnergyBinOriginal = spectrumKeV.photonsPerEnergyBin ;
    spectrumKeV.photonsTotalOriginal =  spectrumKeV.photonsTotal ;
    
    spectrumKeV.photonsPerEnergyBin = computeResultingPhotons( photonsPerEnergyBinKeV, ...
        energyBinLabelsKeV, p.Bowtie.material, mean( spectrumKeV.bowtieThickness(:) ) ) ;
    
    spectrumKeV.photonsTotal =  sum( spectrumKeV.photonsPerEnergyBin ) ;
    spectrumKeV.energyAverage  = sum( spectrumKeV.photonsPerEnergyBin .* spectrumKeV.energyBinLabels ) / spectrumKeV.photonsTotal;
    
    
else
    
    spectrumKeV.useBowtie                   =   false;
    spectrumKeV.photonsTotalOriginal        =  spectrumKeV.photonsTotal ;
    spectrumKeV.photonsPerEnergyBinOriginal = spectrumKeV.photonsPerEnergyBin ;
    
end


%% MV spectrum
[energyBinLabelsMeV, ~, photonsPerEnergyBinMeV, ~, ~, energyAverageMeV] = ...
    readSpectrum([p.Paths.spectraDir p.Spectra.spectrumMeV '.txt'], 'per bin', 1);

photonsPerEnergyBinMeV      = photonsPerEnergyBinMeV * photonsTotalMeV / sum(photonsPerEnergyBinMeV(:));
photonsPerMm2At1mMeV        = photonsTotalMeV * (geomMeV.SDD/1000)^2 / ( geomMeV.detSpacing(1) * geomMeV.detSpacing(end));
photonsPerSteradianMeV      = photonsPerMm2At1mMeV / steradiansPerMm2At1m;

spectrumMeV.energyBinLabels         = energyBinLabelsMeV;
spectrumMeV.photonsPerEnergyBin     = photonsPerEnergyBinMeV;
spectrumMeV.photonsTotal            = photonsTotalMeV ;
spectrumMeV.energyAverage           = energyAverageMeV;
spectrumMeV.photonsPerMm2At1m       = photonsPerMm2At1mMeV ;
spectrumMeV.photonsPerSteradian     = photonsPerSteradianMeV ;
spectrumMeV.DQE                     = p.Detector.detectorConversionEfficiencyMeV;
spectrumMeV.backgroundEvents        = 50;
spectrumMeV.useBowtie               = false;
spectrumMeV.photonsTotalOriginal        =  spectrumMeV.photonsTotal ;
spectrumMeV.photonsPerEnergyBinOriginal = spectrumMeV.photonsPerEnergyBin ;

%% display sprectrum information
if printOut
    
    % print spectra summaries
    pathLengthTissue    = 168;
    pathLengthGold      = 5;
    pathLengthTitanium  = 20;
    pathLengthTissueOnFillingsPath = 146;
    
    photonsEffectiveThroughTissueKeV    = computeResultingPhotons( photonsPerEnergyBinKeV ...
        * p.Detector.detectorConversionEfficiencyKeV, energyBinLabelsKeV, 'Tissue_Soft_ICRU-44', pathLengthTissue);
    
    photonsEffectiveThroughGoldKeV      = computeResultingPhotons( photonsPerEnergyBinKeV ...
        * p.Detector.detectorConversionEfficiencyKeV, energyBinLabelsKeV, 'Tissue_Soft_ICRU-44', pathLengthTissueOnFillingsPath, 'Gold', pathLengthGold );
    
    photonsEffectiveThroughTitaniumKeV  = computeResultingPhotons( photonsPerEnergyBinKeV ...
        * p.Detector.detectorConversionEfficiencyKeV, energyBinLabelsKeV, 'Tissue_Soft_ICRU-44', pathLengthTissueOnFillingsPath, 'Titanium', pathLengthTitanium );
    
    photonsEffectiveThroughTissueMeV    = computeResultingPhotons( photonsPerEnergyBinMeV ...
        * p.Detector.detectorConversionEfficiencyMeV, energyBinLabelsMeV, 'Tissue_Soft_ICRU-44', pathLengthTissue);
    
    photonsEffectiveThroughGoldMeV      = computeResultingPhotons( photonsPerEnergyBinMeV ...
        * p.Detector.detectorConversionEfficiencyMeV, energyBinLabelsMeV, 'Tissue_Soft_ICRU-44', pathLengthTissueOnFillingsPath, 'Gold', pathLengthGold );
    
    photonsEffectiveThroughTitaniumMeV  = computeResultingPhotons( photonsPerEnergyBinMeV ...
        * p.Detector.detectorConversionEfficiencyMeV, energyBinLabelsMeV, 'Tissue_Soft_ICRU-44', pathLengthTissueOnFillingsPath, 'Titanium', pathLengthTitanium );
    
    fprintf('\nSpectra summary:\n');
    fprintf(' Beam | # ph. (P) | # ph. (E) | # ph. (T) | # ph.(Au) | # ph.(Ti) | E_min (keV) | E_max (keV) | E_avg (keV)  \n');
    fprintf('------+-----------+-----------+-----------+-----------+-----------+-------------+-------------+--------------\n');
    fprintf(' keV  | %9.0f | %9.0f | %9.0f | %9.0f | %9.0f |   %7.2f   |   %7.2f   |   %7.2f   \n', ...
        photonsTotalKeV, photonsTotalKeV*p.Detector.detectorConversionEfficiencyKeV, sum(photonsEffectiveThroughTissueKeV), sum(photonsEffectiveThroughGoldKeV), sum(photonsEffectiveThroughTitaniumKeV), energyBinLabelsKeV(1), energyBinLabelsKeV(end), energyAverageKeV);
    
    fprintf(' MeV  | %9.0f | %9.0f | %9.0f | %9.0f | %9.0f |   %7.2f   |   %7.2f   |   %7.2f   \n', ...
        photonsTotalMeV, photonsTotalMeV*p.Detector.detectorConversionEfficiencyMeV, sum(photonsEffectiveThroughTissueMeV), sum(photonsEffectiveThroughGoldMeV), sum(photonsEffectiveThroughTitaniumMeV), energyBinLabelsMeV(1), energyBinLabelsMeV(end), energyAverageMeV);
    fprintf('\n');
    fprintf('# photons (P) are given per pixel at the detector distance for the same current settings.\n');
    fprintf('# photons (E) are given per pixel at the detector distance taking into account the detector efficiency.\n');
    fprintf('# photons (T) are # photons(E) after %gmm of soft tissue.\n', pathLengthTissue);
    fprintf('# photons (Au) are # photons(E) after %gmm of soft tissue and %gmm of Gold fillings each.\n', pathLengthTissueOnFillingsPath, pathLengthGold);
    fprintf('# photons (Ti) are # photons(E) after %gmm of soft tissue and %gmm of Titanium fillings each.\n', pathLengthTissueOnFillingsPath, pathLengthTitanium);
    fprintf('\n');
    fprintf('Additional information for keV imager:\n');
    fprintf('SDD: %.0fmm, Projections: %i, Detector size: %ipx, Pixel size: (%5.3fmm)^2\n', geomKeV.SDD, geomKeV.noViews, p.Geometries.keV.sizeDet, geomKeV.detSpacing);
    %    fprintf('Charge used per image: %4.2fmAs, Total charge: %3.0fmAs, Dose: %5.3fGy\n', mAsPerImage, mAsPerImage*geomKeV.noViews, p.Spectra.doseLimitTotalKeV);
    fprintf('Detector conversion efficiency: %.0f%%\n', 100*p.Detector.detectorConversionEfficiencyKeV);
    fprintf('\n');
    fprintf('Additional information for MeV imager:\n');
    fprintf('SDD: %.0fmm, Projections: %i, Detector size: %ipx, Pixel size: (%5.3fmm)^2\n', geomMeV.SDD, geomMeV.noViews, p.Geometries.MeV.sizeDet, geomMeV.detSpacing);
    %    fprintf('Pulses used per image: %g pulses/image, Total imaging time: %.1fs, Dose: %5.3fGy\n', secondsPerImageMeV*pulseFrequency, secondsPerImageMeV*geomMeV.noViews, p.Spectra.doseLimitTotalMeV);
    fprintf('Detector conversion efficiency: %.0f%%\n', 100*p.Detector.detectorConversionEfficiencyMeV);
    fprintf('\n');
    
end




end