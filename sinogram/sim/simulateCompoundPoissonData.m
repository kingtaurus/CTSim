function intensity = simulateCompoundPossionData( phan, geom, spectrum, sinosDir )
% Compute photon counts data
% Fixed exposure for every projection 
% Detecotor response using Varian high DQE CdWO4 MV detector
%
% input:
%       phan        - phantom parameters
%       geom        - geometry parameters
%       spectrum    - spectrum parameters
%       sinoDire    - sinogram file directory
%       usePolyEnergies ( default true )
% output:
%       intenstiy
%
% Meng Wu at Stanford University
% 2012 - 2013

warning off;

countingNoise           = 1;
spectrumEnergies        = spectrum.energyBinLabels;
spectrumPhotonsPerBin   = spectrum.photonsPerEnergyBinOriginal;
photonsTotal            = spectrum.photonsTotalOriginal;
imgMaterials            = phan.phantomMaterials;

% make sure the caching path exists
if ~exist(sinosDir, 'dir'), mkdir(sinosDir); end
% compute expected number of transmitted photons
fprintf('Compute expected number of transmitted photons ...');

%% check whether cached forward projection exists and get it from file if possible
sinoPhotonFilename = sprintf([sinosDir 'sinoPolyCompoundPoisson_%s_%05.1f_%07.0f.mhd'], ...
    phan.materialMappingName ,spectrum.energyAverage, photonsTotal);
if exist(sinoPhotonFilename, 'file')
    fprintf(' use saved sinogram data. \n');
    intensity = readMetaImage(sinoPhotonFilename);
    return;
end

%% compute expected number of transmitted photons
fprintf('- ray casting ');
tSubStep = tic;

noMaterials = length( phan.materialIndexes );
materialProjection = cell( noMaterials );

for m = 1:noMaterials
    
    progressStr = sprintf('at Material (%i/%i)... ', m, noMaterials );
    
    fprintf(progressStr);
    % check whether cached forward projection exists and get it from file if possible
    sinoMaterialFilename = sprintf([sinosDir 'sinoMaterial%d.mhd'], m );
    
    if exist(sinoMaterialFilename, 'file')
        sinoMaterial = readMetaImage(sinoMaterialFilename);
    else
        % convert material phantom to attenuation image
        imgMaterial = single( imgMaterials == phan.materialIndexes{m} );
        
        if phan.useDensity
            imgMaterial = imgMaterial .* phan.phantomDensities;
        end
        
        % forward projection
        sinoMaterial =  forwardProjectPhantomMex( imgMaterial, phan, geom, 'proj,tf' );
        
        writeMetaImage(sinoMaterial, sinoMaterialFilename);
    end
    materialProjection{m} = sinoMaterial;
end
fprintf('(%is)\n', round(toc(tSubStep)));

%% Bowtie filter
if spectrum.useBowtie 
    sinoBowtie = zeros( size(sinoMaterial), 'single' );
    if length( size(sinoBowtie) ) == 3
        for iv = 1:geom.noViews
            sinoBowtie(:,:,iv) = spectrum.bowtieThickness / 10;
        end
    else
        for iv = 1:geom.noViews
            sinoBowtie(:,iv) = spectrum.bowtieThickness / 10;
        end
    end
else
    sinoBowtie = 0;
end

%% simulate detector noise
intensity           = zeros( size(sinoMaterial), 'single');
effectivePhotons    = zeros( size(sinoMaterial), 'single');

% iterate over all energy bins given and accumulate number of photons
for bin = 1:length(spectrum.energyBinLabels)
    
    photonEnergy = spectrumEnergies(bin);
    
    % print progress figure
    if mod(bin, 10) == 0 || bin == 1
        progressStr = sprintf('at %.1f keV (%i/%i)... ', photonEnergy, bin, length(spectrum.energyBinLabels));
        fprintf(progressStr);
    end
    sinoAtt =  zeros( size(sinoMaterial), 'single');
    
    for m = 1: noMaterials
        materialAtt = materialAttenuation( photonEnergy, phan.materialNames{m});
        sinoAtt = sinoAtt + materialAtt * materialProjection{m} ;
    end
    
    if spectrum.useBowtie
        materialAtt = materialAttenuation( photonEnergy,spectrum.bowtieMaterial);
        sinoAtt = sinoAtt + materialAtt * sinoBowtie ;
    end
    
    % attenuate photons by exponential decay law, take DQE and detectorEnergyRespones into count
    effectivePhotonsBin = spectrumPhotonsPerBin(bin)*exp(-sinoAtt);
    effectivePhotons = effectivePhotons + effectivePhotonsBin;
    
    % simulate detector noise
    if countingNoise
        intensity =  intensity + poissrnd( hDetectorEnergyResponseFunction( effectivePhotonsBin , photonEnergy) ) ...
            * photonEnergy * spectrum.detectorGain;
    else
        intensity = intensity + hDetectorEnergyResponseFunction( effectivePhotonsBin, photonEnergy) ...
            * photonEnergy * spectrum.detectorGain;
    end
    
end

if countingNoise
    intensity = intensity + sqrt( spectrum.electronicNoise ) * randn( size(intensity), 'single' );
end

clear materialProjection sinoAtt;

%% add detector blur
noise = intensity - effectivePhotons;
[effectivePhotons, noise] = addDetectorBlur( effectivePhotons, noise, geom );

%% combine them together
intensity = round( effectivePhotons + noise );
intensity( intensity < 10 ) = 10;
if countingNoise
    writeMetaImage(intensity, sinoPhotonFilename);
end

fprintf('Done\n');

end