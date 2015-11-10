function sinoOut = beamHardeningMaterialCorrection(sinoIn, spectrum, material, maximumThickness )
%
% Copyright (c) 2010-2012 by Meng Wu, Stanford University.
% Modified by Meng Wu 2012.9

if nargin < 4
    maximumThickness = 50;
end

% default material assumed for beam hardening correction is water

if spectrum.energyIntegrating
    photonsPerEnergyBin = spectrum.DQE * spectrum.photonsPerEnergyBin .* spectrum.energyBinLabels * spectrum.detectorGain ;
else
    photonsPerEnergyBin = spectrum.DQE * spectrum.photonsPerEnergyBin * spectrum.detectorGain  ;
end

if spectrum.useBowtie
    photonsPerEnergyBin = computeResultingPhotons( photonsPerEnergyBin, ...
        spectrum.energyBinLabels, spectrum.bowtieMaterial, mean( spectrum.bowtieThickness(:)) );
end

photonsTotal = spectrum.photonsTotal;
energyBinLabels = spectrum.energyBinLabels;

attenuationPoly = materialAttenuation( energyBinLabels, material, photonsPerEnergyBin);
attenuationMono = materialAttenuation( spectrum.energyAverage, material );

% initialize vectors
intensityPolyRef = [];
intensityMonoRef = [];

% compute attenuations for various intersection lengths
% minIntersectionLength = 0.001; % cm
% maxIntersectionLength = 100; % cm
intersectionLengths = logspace(-4, log10(maximumThickness), 512);
for l = intersectionLengths
    intensityPolyBin = sum(photonsPerEnergyBin.*exp(-l*attenuationPoly));
    intensityPolyRef = [intensityPolyRef, intensityPolyBin]; %#ok<AGROW>
    intensityMonoBin = photonsTotal*exp(-l*attenuationMono);
    intensityMonoRef = [intensityMonoRef, intensityMonoBin]; %#ok<AGROW>
end

% compute sinogram values
lineIntegralPoly = [0, -log(intensityPolyRef / photonsTotal)];
lineIntegralMono = [0, -log(intensityMonoRef / photonsTotal)];

if 0
    % plot correction function
    figure('Name', 'Beam Hardening Correction Function', 'NumberTitle', 'off');
    zoom on;
    plot(lineIntegralPoly, lineIntegralMono);
    xlabel('Polychromatic attenuation');
    ylabel(sprintf('Monochromatic attenuation at %.2f keV', energyAverage));
end

sinoIn( sinoIn > max(lineIntegralPoly) - 1e-8 ) = max(lineIntegralPoly) - 1e-8;
sinoIn( sinoIn < 1e-8 ) = 1e-8;

sinoOut = interp1(lineIntegralPoly, lineIntegralMono, sinoIn);
sinoOut(isnan(sinoOut)) = 0;
% zeros really?

end