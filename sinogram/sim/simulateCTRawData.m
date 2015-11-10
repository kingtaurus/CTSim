function [intensity, tubeCurrentProfile] = simulateCTRawData( phan, geom, spectrum, sinosDir, addCountingNoise, sampleNum, upsampleRate )
% Simulate CT raw data
%   most recent function of simulation CT raw data
%   features including:
%       1. Polychomatic x-ray and material based phatom
%       2. focal spot blur, detector noise power spectrum
%       3. Simple possion and compound possion noise
%       4. Automatic exposure control
%
% input:
%       phan        - phantom parameters
%       geom        - geometry parameters
%       spectrum    - spectrum parameters
%       sinoDire    - sinogram file directory
% output:
%       intenstiy
%
% Pervious version:
%   function sinoAtt = simulateAttenuationDataPhantom( phan, geom, spectrum, sinoDir, addNoise )
%   function intensity = simulatePhotonCountingData( phan, geom, spectrum, sinosDir, usePolyEnergies, addCountingNoise )
%   function intensity = simulateCompoundPoissonData( phan, geom, spectrum, sinosDir )
%
% Meng Wu at Stanford University
% 2014


warning off;

if nargin < 5
    addCountingNoise = 1;
end

if nargin < 6
    sampleNum = 1;
end

if nargin < 7
    upsampleRate = 1;
end

%% compute expected number of transmitted photons

fprintf('Compute expected number of transmitted photons ...\n');
% make sure the caching path exists
if ~exist(sinosDir, 'dir'), mkdir(sinosDir); end

% check whether cached forward projection exists and get it from file if possible

sinogramFilename = sprintf([sinosDir 'sinoRawData_%i_%s_kVp%.1f_FLUX%1.2g_BOWT%i_ACE%1.2f_EID%i_NOISE%i.mhd'], ...
    sampleNum, phan.materialMappingName, max(spectrum.energyBinLabels), spectrum.photonsTotalOriginal, ...
    spectrum.useBowtie, spectrum.automaticExposureControl, spectrum.energyIntegrating, addCountingNoise );

tubeCurrentProfileFilename = sprintf([sinosDir 'tubeCurrentProfile_%i_%s_kVp%.1f_FLUX%1.2g_BOWT%i_ACE%1.2f_EID%i_NOISE%i.raw'], ...
    sampleNum, phan.materialMappingName, max(spectrum.energyBinLabels), spectrum.photonsTotalOriginal, ...
    spectrum.useBowtie, spectrum.automaticExposureControl, spectrum.energyIntegrating, addCountingNoise );


if  exist(sinogramFilename, 'file') && exist(tubeCurrentProfileFilename, 'file') %&& 0
    fprintf('\tuse saved sinogram data. \n');
    intensity = readMetaImage(sinogramFilename);
    
    fileID = fopen(tubeCurrentProfileFilename, 'r');
    tubeCurrentProfile = fread(fileID, Inf, 'double');
    fclose(fileID);
    
    return;
end


%% compute expected number of transmitted photons


noMaterials = length( phan.materialIndexes );
materialProjection = cell( noMaterials );


geomProj = geom;
if upsampleRate > 1
    fprintf( '\tUp sampling the detector resolution by %i in forward projection. \n', upsampleRate );
    geomProj.detSize    = geomProj.detSize * upsampleRate;
    geomProj.detOffset  = geomProj.detOffset * upsampleRate;
    geomProj.detSpacing = geomProj.detSpacing / upsampleRate;
end

fprintf('\tRay casting: ');
tSubStep = tic;


for m = 1:noMaterials
    
    progressStr = sprintf('at Material (%i/%i)... ', m, noMaterials );
    
    fprintf(progressStr);
    
    % check whether cached forward projection exists and get it from file if possible
    if upsampleRate == 1
        sinoMaterialFilename = sprintf([sinosDir 'sinoMaterial%d.mhd'], m );
    else
        sinoMaterialFilename = sprintf([sinosDir 'sinoMaterial%dx%i.mhd'], m,  upsampleRate );
    end
    
    if exist(sinoMaterialFilename, 'file') 
        sinoMaterial = readMetaImage(sinoMaterialFilename);
    else
        % convert material phantom to attenuation image
        imgMaterial = single( phan.phantomMaterials == phan.materialIndexes{m} );
        
        if phan.useDensity
            imgMaterial = imgMaterial .* phan.phantomDensities;
        end
        
        % forward projection
        sinoMaterial = forwardProjectPhantomMex( imgMaterial, phan, geomProj, 'proj,tf' );
        
        writeMetaImage(sinoMaterial, sinoMaterialFilename);
        
    end
    
    materialProjection{m} = sinoMaterial;
    
end

fprintf('(%is)\n', round(toc(tSubStep)));

%% Bowtie filter

if spectrum.useBowtie
    
    if upsampleRate == 1
        bowtieThickness = spectrum.bowtieThickness;
    else
        bowtieThickness = loadBowtieFilter( spectrum, geomProj );
    end
    
    sinoBowtie = zeros( size(sinoMaterial), 'single' );
    
    if length( size(sinoBowtie) ) == 3
        for iv = 1:geom.noViews
            sinoBowtie(:,:,iv) = bowtieThickness / 10;
        end
    else
        for iv = 1:geom.noViews
            sinoBowtie(:,iv) = bowtieThickness / 10;
        end
    end
    
else
    sinoBowtie = 0;
end

%% Simulate effective photon count munbers

fprintf('\tSimulation polychromatic sinogram:')
tSubStep = tic;
effectivePhotons = zeros( size(sinoMaterial), 'single');

% iterate over all energy bins given and accumulate number of photons
for bin = 1:length(spectrum.energyBinLabels)
    
    photonEnergy = spectrum.energyBinLabels(bin);
    
    % print progress figure
    if mod(bin, 10) == 0 || bin == 1
        progressStr = sprintf('at %.1f keV (%i/%i)... ', photonEnergy, bin, length(spectrum.energyBinLabels));
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
    effectivePhotonsPerBin = spectrum.DQE * spectrum.photonsPerEnergyBinOriginal(bin) * exp(-sinoAtt);
    
    if spectrum.energyIntegrating
        effectivePhotonsPerBin = effectivePhotonsPerBin * photonEnergy;
    end
    
    if upsampleRate == 1
        effectivePhotons = effectivePhotons + effectivePhotonsPerBin * spectrum.detectorGain;
    else
        if ndims( effectivePhotons ) == 3
            effectivePhotons = effectivePhotons + effectivePhotonsPerBin * spectrum.detectorGain / upsampleRate^2;
        else
            effectivePhotons = effectivePhotons + effectivePhotonsPerBin * spectrum.detectorGain / upsampleRate;
        end
    end
    
end

% downsample and add detector blur
% the order of the following two process are still in question...

effectivePhotons = addDetectorBlur( effectivePhotons, geomProj );

if upsampleRate > 1
    if ndims( effectivePhotons ) == 3
        effectivePhotons = binningDetectorPixels( effectivePhotons, geomProj, [upsampleRate  upsampleRate] );
    else
        effectivePhotons = binningDetectorPixels( effectivePhotons, geomProj, [upsampleRate  1] );
    end
end

clear materialProjection sinoAtt;

fprintf('(%is)\n', round(toc(tSubStep)));


%% Automatic Exposure Control

if spectrum.energyIntegrating
    sourceSpectrum         = spectrum.DQE * spectrum.photonsPerEnergyBin .* spectrum.energyBinLabels * spectrum.detectorGain;
else
    sourceSpectrum         = spectrum.DQE * spectrum.photonsPerEnergyBin * spectrum.detectorGain  ;
end

photonsEffectiveThroughTissue  = sum( computeResultingPhotons( sourceSpectrum, spectrum.energyBinLabels, 'Tissue_Soft_ICRU-44', 250) );

[ effectivePhotons, tubeCurrentProfile ] = automaticExposureControl( effectivePhotons, spectrum, photonsEffectiveThroughTissue, sum( sourceSpectrum ) );

%% simulate detector noise and blur
if spectrum.compoundPoissonNoise
    fprintf('\tWarning: compound Poisson is not currectly implemented.\n')
end

intensity = round( effectivePhotons );

if addCountingNoise
    fprintf('\tSimulating counting noise... ');
    tSubStep = tic;
    
    if ndims( effectivePhotons ) == 3
        
        for iv = 1: size( effectivePhotons, 3 )
            % print progress figure
            if mod(iv, 500) == 0 || iv == 1
                progressStr = sprintf('at (%i/%i) view ... ',  iv, size( effectivePhotons, 3 ) );
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

%% Add Detector Noise Power Spectrum
noise = intensity - effectivePhotons;
noise = addDetectorNoisePowerSpectrum( noise, geom );

intensity = round( effectivePhotons + noise );
intensity( intensity < 1 ) = 1;

%% save for future use
writeMetaImage(intensity, sinogramFilename);

fileID = fopen(tubeCurrentProfileFilename, 'w');
fwrite(fileID, tubeCurrentProfile,'double');
fclose(fileID);

%% Report the Scan dose

maxEnergy = max(spectrum.energyBinLabels( spectrum.photonsPerEnergyBinOriginal > 10 ) );

% use precalculated total number of quanta/mm2/mAs in a spectrum 
if abs( maxEnergy - 140 ) <= 5
    xrnquanta = 3.4762e+06;
elseif abs( maxEnergy - 130 ) <= 5
    xrnquanta = 3.0450e+06;
elseif abs( maxEnergy - 120 ) <= 5
    xrnquanta = 2.6750e+06;
elseif abs( maxEnergy - 110 ) <= 5
    xrnquanta = 2.3243e+06;
elseif abs( maxEnergy - 100 ) <= 5
    xrnquanta = 1.9584e+06;
elseif abs( maxEnergy - 90 ) <= 5
    xrnquanta = 1.5945e+06;
elseif abs( maxEnergy - 80 ) <= 5
    xrnquanta = 1.2412e+06;
elseif abs( maxEnergy - 70 ) < 5
    xrnquanta = 8.7386e+05;
end

fprintf('Photons per pixels: \t%.3g\n', spectrum.photonsTotalOriginal );
fprintf( 'Simulation q / mm2 @1m: \t%.3g \n' , spectrum.photonsPerMm2At1m ); 
fprintf( '1 mAs  q / mm2 @1m: \t\t%.3g \n' , xrnquanta ); 
fprintf( 'Number of projections: \t\t%i \n', length( tubeCurrentProfile ) ); 
nrot = length( tubeCurrentProfile ) / geom.noViews;
fprintf( 'Number of rotation: \t\t%.1f \n', nrot );
fprintf( 'Average mAs per projection: \t%f \n',  mean( tubeCurrentProfile ) * spectrum.photonsPerMm2At1m / xrnquanta );
fprintf( 'Average mAs per rotation: \t%f \n',  sum( tubeCurrentProfile ) * spectrum.photonsPerMm2At1m / xrnquanta / nrot );
fprintf( 'Total mAs: \t\t\t%f \n',  sum( tubeCurrentProfile ) * spectrum.photonsPerMm2At1m / xrnquanta );

fprintf('Done.\n\n');

end