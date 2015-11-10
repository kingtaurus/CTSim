function [CNR photons] = spectrumCNR( spectrum, thickness, ltissue, lbone, lmetal  )


% Compute attenuations lookup table for two sprectrums
materialsDir = 'physicsdata/materials/';

[E_tissue1, mu_tissue1] = readMaterialAttenuation('Adipose_Tissue_ICRU-44', materialsDir);

[E_tissue2, mu_tissue2] = readMaterialAttenuation('Tissue_Soft_ICRU-44', materialsDir);

[E_bone, mu_bone] = readMaterialAttenuation('Bone_Cortical_ICRU-44', materialsDir);

[E_metal, mu_metal] = readMaterialAttenuation('Gold', materialsDir);


muTissue1 = spectralAttenuations(spectrum.energyBinLabels, E_tissue1, mu_tissue1);

muTissue2 = spectralAttenuations(spectrum.energyBinLabels, E_tissue2, mu_tissue2);

muBone = spectralAttenuations(spectrum.energyBinLabels, E_bone, mu_bone);

muMetal = spectralAttenuations(spectrum.energyBinLabels, E_metal, mu_metal);

intensity = hDetectorEnergyResponseFunction(spectrum.photonsPerEnergyBin ...
    .*exp(- ltissue * muTissue2 - lbone * muBone - lmetal * muMetal ), spectrum.photonsPerEnergyBin);
	

constrast = exp( thickness * ( muTissue2 - muTissue1 ) ) .*  intensity;

CNR = log( sum ( constrast ) ) ./ sqrt(  1 / sum( intensity ) )
%CNR = sum( cnrs .*  intensity / sum(intensity));

photons = sum( intensity )



figure
plot( spectrum.energyBinLabels, spectrum.photonsPerEnergyBin / sum( spectrum.photonsPerEnergyBin )); hold on;
plot( spectrum.energyBinLabels, intensity / sum( intensity ), 'r');


figure
plot(spectrum.energyBinLabels, muTissue1, '-' ); hold on;
plot(spectrum.energyBinLabels, muTissue2, '-.' ); hold on;
plot(spectrum.energyBinLabels, muBone, '--' );

% 
% figure
% loglog(spectrum.energyBinLabels, muMetal, '--' );
% grid on;
% axis tight;



end