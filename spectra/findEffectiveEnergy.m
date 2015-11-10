function energyEff = findEffectiveEnergy( energyBinLabels, photonsPerEnergyBin, thickness )


muWater = materialAttenuation( energyBinLabels, 'water' );

photonsResulting = computeResultingPhotons(photonsPerEnergyBin, energyBinLabels, 'water', thickness * 10 );

attenApprox = log( sum(photonsPerEnergyBin) / sum( photonsResulting ) ) / thickness ; 

energyEff = interp1(  muWater, energyBinLabels,attenApprox);

end