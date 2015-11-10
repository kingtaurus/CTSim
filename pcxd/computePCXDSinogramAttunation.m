function sinoAtt = computePCXDSinogramAttunation( sinoPhotonCounts, energyBinSizes )
% Compute sinogram attenuation for PCDX by -log(I/I_0)
% input:
%       intensity - photon counts data of sinogram
%       spectrm - x-ray spectrum
% output:
%       sinoAtt - sinogram attenuation
%
% Modified by Meng Wu 2014.3

% normalize by expected intensity for unattenuated source spectrum photon densities

fprintf('Computing atteunation for PCXD ...');

dummyLowPhotonCount = 10;
noDetectorEnergyBins = length( energyBinSizes );
sinoAtt = cell( noDetectorEnergyBins, 1 );

% convert back to attenuation integrals (preventing a NaN by ensuring that at least "dummyLowPhotonCount photons are detected" at any pixel)


for ib = 1 : noDetectorEnergyBins
    intensity = sinoPhotonCounts{ib};
    intensity = max( intensity,  dummyLowPhotonCount);
    sinoAtt{ib} = -log(intensity / energyBinSizes(ib) );
    
end

fprintf('Done\n\n');


end