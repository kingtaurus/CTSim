function [ Ybar, nablaYsoft, nablaYbone, tableSoft, tableBone, muEffSoft, muEffBone, muEffMetal] ...
    = computePolychromaticAttLookupTable( spectrum, material1, material2, material3 )

if nargin < 2
    material1 = 'Tissue_Soft_ICRU-44';
    material2 = 'Bone_Cortical_ICRU-44';
end

if nargin < 4
    material3 = 'Fe';
end

tableSize = 512;
minLengthSoft = -1e-8;
maxLengthSoft = 60;
minLengthBone = -1e-8;
maxLengthBone = 20;

spectrumEnergies        = spectrum.energyBinLabels;
spectrumPhotonsPerBin   = spectrum.photonsPerEnergyBin / sum( spectrum.photonsPerEnergyBin ) ...
    * spectrum.DQE * sum( spectrum.photonsPerEnergyBinOriginal ) ;

[ muPoly_Soft, muEffSoft ] = materialAttenuation( spectrumEnergies, material1, spectrumPhotonsPerBin );
[ muPoly_Bone, muEffBone ] = materialAttenuation( spectrumEnergies, material2, spectrumPhotonsPerBin );
[ ~, muEffMetal ] = materialAttenuation( spectrumEnergies, material3, spectrumPhotonsPerBin );

lenSoft = linspace( minLengthSoft, maxLengthSoft, tableSize );
lenBone = linspace( minLengthBone, maxLengthBone, tableSize );

[ tableSoft, tableBone ] = meshgrid( lenSoft, lenBone );

Ybar = zeros( tableSize );
nablaYsoft = zeros( tableSize );
nablaYbone = zeros( tableSize );

%compute the lookup table for Y bar
for i = 1 : tableSize
    for j = 1 : tableSize
        
        ls = tableSoft(i, j);
        
        lb = tableBone(i, j);
        
        Ybar(i,j) =  sum( spectrumPhotonsPerBin.*exp(- ls * muPoly_Soft - lb * muPoly_Bone ));
        
        nablaYsoft(i,j) =  - sum( muPoly_Soft.* spectrumPhotonsPerBin.*exp(- ls * muPoly_Soft - lb * muPoly_Bone ));
        
        nablaYbone(i,j) =  - sum( muPoly_Bone.* spectrumPhotonsPerBin.*exp(- ls * muPoly_Soft - lb * muPoly_Bone ));
        
    end
end

end