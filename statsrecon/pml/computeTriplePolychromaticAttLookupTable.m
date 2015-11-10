function [ Ybar, nablaYsoft, nablaYbone, nablaYMetal, tableSoft, tableBone, tableMetal, ...
    muEff_Soft, muEff_Bone, muEff_Metal] = computeTriplePolychromaticAttLookupTable( spectrum , metalMaterial, material1, material2  )

if nargin < 4
    material1 = 'Tissue_Soft_ICRU-44';
    material2 = 'Bone_Cortical_ICRU-44';
end

tableSize = 256;
minLengthSoft = 0;
maxLengthSoft = 60;
minLengthBone = 0;
maxLengthBone = 30;
minLengthMetal = -0.1;
maxLengthMetal = 10;

spectrumEnergies        = spectrum.energyBinLabels;
%spectrumPhotonsPerBin   = spectrum.photonsPerEnergyBin * spectrum.DQE;
spectrumPhotonsPerBin   = hDetectorEnergyResponseFunction( spectrum.photonsPerEnergyBin , spectrumEnergies) .* spectrumEnergies / mean(spectrumEnergies) ;

[ muPoly_Soft, muEff_Soft ] = materialAttenuation( spectrumEnergies, material1, spectrumPhotonsPerBin );
[ muPoly_Bone, muEff_Bone ] = materialAttenuation( spectrumEnergies, material2, spectrumPhotonsPerBin );
[ muPoly_Metal, muEff_Metal ] = materialAttenuation( spectrumEnergies, metalMaterial, spectrumPhotonsPerBin );

lenSoft = linspace( minLengthSoft, maxLengthSoft, tableSize );
lenBone = linspace( minLengthBone, maxLengthBone, tableSize );
lenMetal = linspace( minLengthMetal, maxLengthMetal, tableSize/2 );

[ tableSoft, tableBone, tableMetal ] = meshgrid( lenSoft, lenBone, lenMetal );

Ybar = zeros( size(tableSoft) );
nablaYsoft = zeros( size(tableSoft) );
nablaYbone = zeros( size(tableSoft) );
nablaYMetal = zeros( size(tableSoft) );


%compute the lookup table for Y bar
for i = 1 : size(tableSoft, 1)
    for j = 1 : size(tableSoft, 2)
        for k = 1 : size(tableSoft, 3)
            
            ls = tableSoft(i, j, k);
            
            lb = tableBone(i, j, k);
            
            lm = tableMetal(i, j, k);
            
            
            Ybar(i,j, k) =  sum( spectrumPhotonsPerBin.*exp(- ls * muPoly_Soft - lb * muPoly_Bone - lm * muPoly_Metal ));
            
            nablaYsoft(i,j, k) =  - sum(  muPoly_Soft .* spectrumPhotonsPerBin.*exp(- ls * muPoly_Soft - lb * muPoly_Bone - lm * muPoly_Metal ));
            
            nablaYbone(i,j, k) =  - sum(  muPoly_Bone .* spectrumPhotonsPerBin.*exp(- ls * muPoly_Soft - lb * muPoly_Bone - lm * muPoly_Metal ));
            
            nablaYMetal(i,j, k) = - sum(  muPoly_Metal .* spectrumPhotonsPerBin.*exp(- ls * muPoly_Soft - lb * muPoly_Bone - lm * muPoly_Metal ));
            
            
            
        end
    end
end


end