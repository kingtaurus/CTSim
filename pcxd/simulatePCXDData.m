function [ sinoPhotonCounts, energyBinSizes ] = simulatePCXDData( phan, geom, spectrum, sinosDir, binThresholds, countingNoise )
% simulate photon counts data
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


if nargin < 6
    countingNoise = true;
end

warning off;

spectrumEnergies        = spectrum.energyBinLabels;
spectrumPhotonsPerBin   = spectrum.photonsPerEnergyBinOriginal;
photonsTotal            = spectrum.photonsTotalOriginal;
imgMaterials            = phan.phantomMaterials;
backgroundEvents        = spectrum.backgroundEvents;

materialIndexes         = phan.materialIndexes;
materialNames           = phan.materialNames;

noEnergyBins = length(spectrum.energyBinLabels);

% make sure the caching path exists
if ~exist(sinosDir, 'dir'), mkdir(sinosDir); end
% compute expected number of transmitted photons
fprintf([ 'Simulate PCXD data with bin thresholds at ' num2str(binThresholds) ' keV ... \n'] );

%% compute expected number of transmitted photons
fprintf('\t- ray casting ');
tSubStep = tic;

noMaterials = length( materialIndexes );
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
        imgMaterial = single( imgMaterials == materialIndexes{m} );
        
        if phan.useDensity
            imgMaterial = imgMaterial .* phan.phantomDensities;
        end
        
        % forward projection
        sinoMaterial =  forwardProjectPhantomMex( imgMaterial, phan, geom, 'proj,rd' );
        
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
fprintf('\t- compute effective photon numbers ')
    
noDetectorEnergyBins = length( binThresholds );
sinoPhotonCounts = cell( noDetectorEnergyBins, 1 );
energyBinSizes = zeros( noDetectorEnergyBins, 1 );
effectivePhotons = zeros( size(sinoMaterial), 'single');

ibp = 1;
% iterate over all energy bins given and accumulate number of photons
for bin = 1:noEnergyBins
    

    ib = find( binThresholds < spectrumEnergies(bin), 1, 'last' );
    if ib ~= ibp
        sinoPhotonCounts{ibp} = effectivePhotons;
        effectivePhotons = zeros( size(sinoMaterial), 'single');
        ibp = ib;
    end
    
    % print progress figure
    if mod(bin, 10) == 0 || bin == 1
        progressStr = sprintf('at %.1f keV (%i/%i)... ', spectrumEnergies(bin), bin, noEnergyBins);
        fprintf(progressStr);
    end
    
    sinoAtt =  zeros( size(sinoMaterial), 'single');
    
    for m = 1: noMaterials
        materialAtt = materialAttenuation( spectrumEnergies(bin), materialNames{m});
        sinoAtt = sinoAtt + materialAtt * materialProjection{m} ;
    end
    
    if spectrum.useBowtie
        materialAtt = materialAttenuation( spectrumEnergies(bin),spectrum.bowtieMaterial);
        sinoAtt = sinoAtt + materialAtt * sinoBowtie ;
    end
    
    % attenuate photons by exponential decay law, take DQE and detectorEnergyRespones into count
    effectivePhotons = effectivePhotons + spectrum.DQE * spectrumPhotonsPerBin(bin)*exp(-sinoAtt);
   
    energyBinSizes(ib) = energyBinSizes(ib) + spectrum.DQE * spectrumPhotonsPerBin(bin);
    
end
%last bin
sinoPhotonCounts{ib} = effectivePhotons;

clear materialProjection sinoAtt;

fprintf('(%is)\n', round(toc(tSubStep)));

%% simulate detector noise
if countingNoise
    fprintf('\t- simulating counting noise... ');
    tSubStep = tic;
    
    for ib = 1: noDetectorEnergyBins
        
        fprintf('at bin (%i/%i)... ',  ib, noDetectorEnergyBins);
        
        effectivePhotons = sinoPhotonCounts{ib};
        if length( size(effectivePhotons) ) == 2
            effectivePhotons =  poissrnd( effectivePhotons ) ;
        else
            for iv = 1:size(effectivePhotons, 3)
                effectivePhotons(:,:,iv) = poissrnd( effectivePhotons(:,:,iv) );
            end
        end
        sinoPhotonCounts{ib} = effectivePhotons;
    end

    fprintf('(%is)\n', round(toc(tSubStep)));
else
    fprintf('- NOT simulating counting noise\n');
end

fprintf('Done!\n\n');

end