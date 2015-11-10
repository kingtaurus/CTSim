function [muSoft, muBone] = getSoftBoneAttenuation( spectrum )

[~, muSoft] = materialAttenuation( spectrum.energyBinLabels , 'Tissue_Soft_ICRU-44', spectrum.photonsPerEnergyBin);
[~, muBone] = materialAttenuation( spectrum.energyBinLabels , 'Bone_Cortical_ICRU-44', spectrum.photonsPerEnergyBin);

end