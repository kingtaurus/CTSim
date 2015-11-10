function spectrum= loadSpectraCT(p, geom, numPhotonsPerPixel )
% Load X-ray spectra
%   inputs:
%       p - configuration paramters
%       geom - geometry paramters
%       numPhotonsPerPixel 
%   output:
%       spectrum - spectrum parameter
%
% Copyright (c) 2010-2012 by Andreas Keil, Stanford University.
% Modified by Meng Wu at 2012.9

if nargin < 3
    numPhotonsPerPixel = 1e6;
end

% read spectrum from the file
[energyBinLabels, ~, photonsPerEnergyBin ]  = readSpectrum([p.Paths.spectraDir p.Spectra.spectrum '.txt'], 'per bin', 1);
spectrum.energyBinLabels                    = energyBinLabels;

%% scale the spectrum if using the energy integrating detector
spectrum.energyIntegrating      = p.Detector.energyIntegrating;
spectrum.compoundPoissonNoise   = p.Detector.compoundPoissonNoise;
if spectrum.energyIntegrating
    spectrum.detectorGain   = sum( photonsPerEnergyBin ) / ( sum(energyBinLabels .* photonsPerEnergyBin ) );
else
    spectrum.detectorGain	= 1;
end

photonsPerEnergyBin             = photonsPerEnergyBin * numPhotonsPerPixel / sum(photonsPerEnergyBin(:));
spectrum.photonsPerEnergyBin    = photonsPerEnergyBin;
spectrum.photonsTotal           = numPhotonsPerPixel;

%% other spectrum parameters
if spectrum.energyIntegrating
    spectrum.energyAverage = findEffectiveEnergy( energyBinLabels, photonsPerEnergyBin .* energyBinLabels , 30 ); 
else
    spectrum.energyAverage = findEffectiveEnergy( energyBinLabels, photonsPerEnergyBin , 30 ); 
end

spectrum.DQE                    = p.Detector.detectorConversionEfficiency;
spectrum.electronicNoise        = 50;
spectrum.useBowtie              = false;
spectrum.automaticExposureControl   = p.Spectra.automaticExposureControl;


%% Compute standard x-ray intensity measure

steradiansPerMm2At1m            = solidAnglePyramidFromDistanceAndBaseLengths(1, 0.001, 0.001);
perPixel2PerMm2At1m             = (geom.SDD/1000)^2 / ( geom.detSpacing(1) * geom.detSpacing(end));
photonsPerMm2At1m               = numPhotonsPerPixel * perPixel2PerMm2At1m;
photonsPerSteradian             = photonsPerMm2At1m / steradiansPerMm2At1m;
spectrum.maxPhotonsPerPixel     = p.Spectra.maximumIntensity / perPixel2PerMm2At1m;



spectrum.photonsPerMm2At1m      = photonsPerMm2At1m;
spectrum.photonsPerSteradian    = photonsPerSteradian;

%% add bowtie filter into kV spectrum

spectrum.Bowtie = p.Bowtie;

if ~strcmpi(p.Bowtie.shapeType, 'none')
    
    spectrum.useBowtie       = true;
    spectrum.bowtieThickness = loadBowtieFilter( p, geom );
    spectrum.bowtieMaterial  =  p.Bowtie.material;
    spectrum.flatFieldRatio  = zeros( size(spectrum.bowtieThickness) );
    
    for i = 1:size( spectrum.flatFieldRatio, 2)
        
        if spectrum.energyIntegrating
            
            spectrum.flatFieldRatio(:,i) = sum( energyBinLabels.* computeResultingPhotons( photonsPerEnergyBin, ...
                energyBinLabels, p.Bowtie.material, spectrum.bowtieThickness(1,i) ) ) / sum( energyBinLabels .* photonsPerEnergyBin);

        else
            
            spectrum.flatFieldRatio(:,i) = sum( computeResultingPhotons( photonsPerEnergyBin, ...
                energyBinLabels, p.Bowtie.material, spectrum.bowtieThickness(1,i) ) ) / numPhotonsPerPixel;
        end
        
    end
    
    spectrum.photonsPerEnergyBinOriginal    = spectrum.photonsPerEnergyBin ;
    spectrum.photonsTotalOriginal           =  spectrum.photonsTotal ;
    
    spectrum.photonsPerEnergyBin = computeResultingPhotons( photonsPerEnergyBin, ...
        energyBinLabels, p.Bowtie.material, min(spectrum.bowtieThickness(:))  ) ;
    
    spectrum.photonsTotal 	= sum( spectrum.photonsPerEnergyBin ) ;
    
    % The average energy got shifted too. but i am not compensate it, since
    % it wasn't right at the beginning.
else
    
    spectrum.useBowtie                   = false;
    spectrum.photonsTotalOriginal        = spectrum.photonsTotal ;
    spectrum.photonsPerEnergyBinOriginal = spectrum.photonsPerEnergyBin ;
    
end


%% print out summary
printOutSpectrumInfo( spectrum )

fprintf('Done.\n\n');



end