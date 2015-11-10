function intensity = simulateSimplePoissonData( phan, geom, spectrum, sinosDir, usePolyEnergies, addCountingNoise )
% simulate photon counts data
% fixed tube mAs
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

if nargin < 5
    usePolyEnergies = true;
end


if nargin < 6
    addCountingNoise = true;
end


warning off;

spectrumEnergies        = spectrum.energyBinLabels;
spectrumPhotonsPerBin   = spectrum.photonsPerEnergyBinOriginal;
photonsTotal            = spectrum.photonsTotalOriginal;
imgMaterials            = phan.phantomMaterials;

noEnergyBins = length(spectrum.energyBinLabels);

% make sure the caching path exists
if ~exist(sinosDir, 'dir'), mkdir(sinosDir); end
% compute expected number of transmitted photons
fprintf('Compute expected number of transmitted photons ...');


%% check whether cached forward projection exists and get it from file if possible
if usePolyEnergies
    sinoPhotonFilename = sprintf([sinosDir 'sinoPoly_%s_kVp%3.0f_FLUX%1.2g_BOWT%i_ACE%i_EID%i_NOISE%i.mhd'], ...
        phan.materialMappingName, max(spectrum.energyBinLabels), spectrum.photonsTotal, ...
        spectrum.useBowtie, spectrum.automaticExposureControl, spectrum.energyIntegrating, addCountingNoise );
else
    sinoPhotonFilename = sprintf([sinosDir 'sinMono_%s_kVp%3.0f_FLUX%1.2g_BOWT%i_ACE%i_EID%i_NOISE%i.mhd'], ...
        phan.materialMappingName, max(spectrum.energyBinLabels), spectrum.photonsTotal, ...
        spectrum.useBowtie, spectrum.automaticExposureControl, spectrum.energyIntegrating, addCountingNoise );
end

if  exist(sinoPhotonFilename, 'file') && 0
    fprintf(' use saved sinogram data. \n');
    intensity = readMetaImage(sinoPhotonFilename);
    return;
end


%% compute expected number of transmitted photons
fprintf('- ray casting ');
tSubStep = tic;

noMaterials = length(  phan.materialIndexes );
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

%%
effectivePhotons = zeros( size(sinoMaterial), 'single');

if usePolyEnergies
    
    % iterate over all energy bins given and accumulate number of photons
    for bin = 1:noEnergyBins
        
        photonEnergy = spectrumEnergies(bin);
        
        % print progress figure
        if mod(bin, 10) == 0 || bin == 1
            progressStr = sprintf('at %.1f keV (%i/%i)... ', photonEnergy, bin, noEnergyBins);
            fprintf(progressStr);
        end
        sinoAtt =  zeros( size(sinoMaterial), 'single');
        
        for m = 1: noMaterials
            materialAtt = materialAttenuation( photonEnergy,  phan.materialNames{m});
            sinoAtt = sinoAtt + materialAtt * materialProjection{m} ;
        end
        
        if spectrum.useBowtie
            materialAtt = materialAttenuation( photonEnergy, spectrum.bowtieMaterial);
            sinoAtt = sinoAtt + materialAtt * sinoBowtie ;
        end
        
        % attenuate photons by exponential decay law, take DQE and detectorEnergyRespones into count
        effectivePhotonsPerBin = spectrum.DQE * spectrumPhotonsPerBin(bin) * exp(-sinoAtt);
        
        if spectrum.energyIntegrating
            effectivePhotonsPerBin = effectivePhotonsPerBin * photonEnergy;
        end
        
        effectivePhotons = effectivePhotons + effectivePhotonsPerBin * spectrum.detectorGain;
    end
    
else
    
    sinoAtt =  zeros( size(sinoMaterial), 'single');
    
    for m = 1: noMaterials
        [ ~, materialAtt] = materialAttenuation( spectrumEnergies, phan.materialNames{m}, spectrumPhotonsPerBin);
        sinoAtt = sinoAtt + materialAtt * materialProjection{m} ;
    end
    
    if spectrum.useBowtie
        [~, materialAtt] = materialAttenuation( spectrumEnergies,spectrum.bowtieMaterial, spectrumPhotonsPerBin);
        sinoAtt = sinoAtt + materialAtt * sinoBowtie ;
    end
    
    % attenuate photons by exponential decay law, take DQE and detectorEnergyRespones into count
    effectivePhotons = spectrum.DQE * photonsTotal * exp( -sinoAtt );
    
end

clear materialProjection sinoAtt;

fprintf('(%is)\n', round(toc(tSubStep)));

%% simulate detector noise

intensity = round( effectivePhotons );

if addCountingNoise
    fprintf('\tSimulating counting noise... ');
    tSubStep = tic;
    
    if ndims( effectivePhotons ) == 3
        
        for iv = 1: size( effectivePhotons, 3 )    
            % print progress figure
            if mod(iv, 500) == 0 || iv == 1
                progressStr = sprintf('at %i/%i view ... ',  iv, size( effectivePhotons, 3 ) );
                fprintf(progressStr);
            end

            view = effectivePhotons(:,:,iv);
            intensity(:,:,iv) =  poissrnd( view ) + sqrt( spectrum.electronicNoise ) * randn( size(view), 'single' ); 
        end
    else
        intensity =  poissrnd( effectivePhotons ) + sqrt( spectrum.electronicNoise ) * randn( size(effectivePhotons), 'single' );
    end
    
    fprintf('(%is)\n', round(toc(tSubStep)));
else
    fprintf('\tNot Simulating counting noise\n');
    
end

noise = intensity - effectivePhotons;

[effectivePhotons, noise] = addDetectorBlur( effectivePhotons, noise, geom );

%% combine them together

intensity = round( effectivePhotons + noise );
intensity( intensity < 10 ) = 10;
if addCountingNoise
    writeMetaImage(intensity, sinoPhotonFilename);
end

end