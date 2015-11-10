function sinoAtt = computeSimplePoissonSinogramAttunation( intensity, spectrum )
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

dummyLowPhotonCount = 10;

intensity = max(  intensity,  dummyLowPhotonCount);

if length( spectrum.energyBinLabels ) == 1
    
    fprintf('\tmonochromatic spectrum with photon counting detector... \n');
    spectrumPhotons     = spectrum.photonsTotalOriginal;
    sourceIntensity     = spectrum.DQE * spectrumPhotons;
    
else
    fprintf('\tploychromatic spectrum ');
    spectrumPhotonsPerBin   = spectrum.photonsPerEnergyBinOriginal;
    
    if spectrum.energyIntegrating
        fprintf( 'with energy integrating detector ... \n')
        sourceIntensity         = spectrum.DQE * sum(spectrumPhotonsPerBin .* spectrum.energyBinLabels ) * spectrum.detectorGain ;
    else
        fprintf( 'with photon counting detector ... \n')
        sourceIntensity         = spectrum.DQE * sum(spectrumPhotonsPerBin) * spectrum.detectorGain  ;
    end
end

intensityRatio          = intensity / sourceIntensity;

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

% log
sinoAtt = -log(intensityRatio);

end