function imgOut = attMono2Poly( imgIn, spectrum )

% Compute attenuations lookup table for two sprectrums
materialsDir = 'physicsdata/materials/';

materials = {'Air_Dry_near_sea_level', ...
    'Adipose_Tissue_ICRU-44', ...
    'Water_liquid', ...
    'Muscle_Skeletal_ICRU-44', ...
    'Tissue_Soft_ICRU-44', ...
    'B-100_Bone-Equivalent_Plastic', ...
    'Bone_Cortical_ICRU-44'};

[E, mu] = readMaterialAttenuations(materials, materialsDir);

noMaterials = length(materials);

muPoly = zeros(noMaterials, 1);
muMono = zeros(noMaterials, 1);

for m = 1 : noMaterials
    
    muMono(m) = spectralAttenuation(spectrum.energyAverage, 1, E{m}, mu{m});
    
    muPoly(m) = spectralAttenuation(spectrum.energyBinLabels, ...
        spectrum.photonsPerEnergyBin, E{m}, mu{m});
    
end

% make sure no pixel is outside the lookup table
imgIn( imgIn < muMono(2) * 0.5 ) = muMono(1);
imgIn( imgIn > muMono(end) ) =  muMono(end);

imgOut = interp1( muMono, muPoly, imgIn);



end