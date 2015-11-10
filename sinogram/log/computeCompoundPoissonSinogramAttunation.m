function sinoAtt = computeCompoundPoissonSinogramAttunation( intensity, spectrum )
% Compute sinogram attenuation integrals by -log(I/I_0)
% input:
%       intensity - photon counts data of sinogram
%       spectrm - x-ray spectrum
% output:
%       sinoAtt - sinogram attenuation
%
% Copyright (c) 2010-2012 by Andreas Keil, Stanford University.
% Modified by Meng Wu 2012.9

% normalize by expected intensity for unattenuated source spectrum photon densities

dummyLowPhotonCount = 100;

spectrumEnergies        = spectrum.energyBinLabels;
spectrumPhotonsPerBin   = spectrum.photonsPerEnergyBinOriginal;

sourceIntensity         = sum( hDetectorEnergyResponseFunction(spectrumPhotonsPerBin, spectrumEnergies) .* spectrumEnergies * spectrum.detectorGain );
intensityRatio          = intensity / sourceIntensity;

fprintf('\tenergy integrating detector with compound poisson model ...\n')
if spectrum.useBowtie
    if length( size(intensityRatio) ) == 3
        for iv = 1: size(intensityRatio, 3 )
            intensityRatio(:,:,iv) = intensityRatio(:,:,iv) ./ spectrum.flatFieldRatio;
        end
    else
        spectrum.flatFieldRatio = spectrum.flatFieldRatio';
        for iv = 1:size(intensityRatio, 2 )
            intensityRatio(:,iv) = intensityRatio(:,iv) ./ spectrum.flatFieldRatio;
        end
    end
end

% convert back to attenuation integrals (preventing a NaN by ensuring that at least "dummyLowPhotonCount photons are detected" at any pixel)

dummyIntensityRatio = dummyLowPhotonCount / sourceIntensity;
sinoAtt             = - log(max(intensityRatio, dummyIntensityRatio*ones(size(intensityRatio))));

end