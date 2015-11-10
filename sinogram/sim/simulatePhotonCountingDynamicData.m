function intensity = simulatePhotonCountDynamicSinogram( phan, geom, spectrum, sinosDir, TempContrast )
% function sinoAtt = computePloyAttenuationSinogram(...)
%   computes the noisy attenuation intergral sinograms by transmitting
%   photons through a voxelized materials phantom, taking into account
%   noise sources, and relating it to input intensities again (-ln(I/I_0)).
%
%
% Copyright (c) 2010-2012 by Andreas Keil, Stanford University.
% Modified by Meng Wu 2012.9

spectrumEnergies        = spectrum.energyBinLabels;
spectrumPhotonsPerBin   = spectrum.photonsPerEnergyBin;
imgMaterials            = phan.phantomMaterials;
backgroundEvents        = spectrum.backgroundEvents;

materialIndexes         = phan.materialIndexes;
materialNames           = phan.materialNames;
materialsDir            = phan.materialsDir;


noViews                 = geom.noViews;
sizeDet                 = geom.detSize;
spacingDet              = geom.detSpacing;
countingNoise           = false;

noEnergyBins = length(spectrum.energyBinLabels);


% make sure the caching path exists
if ~exist(sinosDir, 'dir'), mkdir(sinosDir); end

% compute expected number of transmitted photons
fprintf('Compute expected number of transmitted photons ...');

% compute expected number of transmitted photons
fprintf('- ray casting ');
tSubStep = tic;


noMaterials = length( materialIndexes );
materialProjection = zeros(sizeDet, noViews, noMaterials );

for m = 1:noMaterials
    
    progressStr = sprintf('at Material (%i/%i)... ', m, noMaterials );
    
    fprintf(progressStr);
    % check whether cached forward projection exists and get it from file if possible
    sinoMaterialFilename = sprintf([sinosDir 'sinoPolyMaterial%d.mhd'], m );
    
    if exist(sinoMaterialFilename, 'file')
        sinoMaterial = readMetaImage(sinoMaterialFilename);
    else
        % convert material phantom to attenuation image
        imgMaterial = single( imgMaterials == materialIndexes{m} );
        
       % forward projection
        sinoMaterial =  forwardProjectPhantomMex( imgMaterial, phan, geom, 'proj,tf' );
        
        writeMetaImage(sinoMaterial, sinoMaterialFilename );
    end
    
    
    materialProjection(:,:,m) = sinoMaterial;
    
end
fprintf('(%is)\n', round(toc(tSubStep)));

% init sinogram
sinoPhotons = zeros(sizeDet, noViews, noEnergyBins);
% iterate over all energy bins given and accumulate number of photons
for bin = 1:noEnergyBins
    % print progress figure
    progressStr = sprintf('at %.1f keV (%i/%i)... ', spectrumEnergies(bin), bin, noEnergyBins);
    fprintf(progressStr);
    sinoAtt = zeros( sizeDet, noViews );
    
    for m = 1: noMaterials     
        if m < noMaterials
            materialAtt = materialAttenuation( spectrumEnergies(bin), materialNames{m});
            
            sinoAtt = sinoAtt + materialAtt * squeeze( materialProjection(:,:,m) ) ;
            
        else
            materialAtt = materialAttenuation( spectrumEnergies(bin), materialNames{2});
            
            contrastAtt = materialAttenuation( spectrumEnergies(bin), materialNames{m});
            
            sinoAtt = sinoAtt + materialAtt * squeeze( materialProjection(:,:,m) ) ;
            
            for i = 1 : geom.noViews
            
            sinoAtt(:,i) = sinoAtt(:,i) + TempContrast(i) *  contrastAtt * squeeze( materialProjection(:,i,m) ) ;
            
            end
        end
        
    end
    
    % attenuate photons by exponential decay law
    sinoPhotons(:, :, bin) = spectrumPhotonsPerBin(bin)*exp(-sinoAtt);
    
end
fprintf('(%is)\n', round(toc(tSubStep)));

% take DQE(0) into account
fprintf('- taking detector conversion efficiency (DQE(0)) into account... ');
effectivePhotons = geom.DQE * sinoPhotons;
clear sinoPhotons;
fprintf('\n');

% simulate detector noise
if countingNoise
    fprintf('- simulating counting noise... ');
    tSubStep = tic;
    
    detectedPhotons = poissrnd( effectivePhotons );
    
    fprintf('(%is)\n', round(toc(tSubStep)));
else
    fprintf('- NOT simulating counting noise\n');
    detectedPhotons = effectivePhotons;
    clear effectivePhotons;
end

% take energy-dependant detector response into account
fprintf('- simulating energy-dependant response... ');
intensity = zeros(size(detectedPhotons));
% compute intensity per energy bin
for b = 1:noEnergyBins
    intensity(:, :, b) = hDetectorEnergyResponseFunction(squeeze(detectedPhotons(:, :, b)), spectrumEnergies(b));
end
% sum up responses over all energy bins
intensity = squeeze(sum(intensity, 3)) + backgroundEvents;
% clean up
clear detectedPhotons;
fprintf('\n');

end