function crlb = cramerRaoLowerBoundPhotonCountingDetector( spectrum, Ttissue, Tbone, Tmetal  )

energyBinLabels = spectrum.energyBinLabels;
photonsPerEnergyBin = spectrum.photonsPerEnergyBin;

photonsPerEnergyBin = photonsPerEnergyBin * 1e4 / sum( photonsPerEnergyBin );

photonsPerEnergyBin = hDetectorAbsorbResponseFunction( photonsPerEnergyBin,  energyBinLabels);

muSoftPoly = materialAttenuation( energyBinLabels, 'Tissue_Soft_ICRU-44', photonsPerEnergyBin );
muBonePoly = materialAttenuation( energyBinLabels, 'Bone_Cortical_ICRU-44', photonsPerEnergyBin );
muMetalPoly  = materialAttenuation( energyBinLabels, 'Gold', photonsPerEnergyBin );


q = muSoftPoly * Ttissue + muBonePoly * Tbone + muMetalPoly * Tmetal;


lambda = sum( photonsPerEnergyBin .* exp( - q ) );

dlambda = sum( muSoftPoly .* photonsPerEnergyBin .* exp( - q ) );

crlb = lambda / dlambda^2 ;