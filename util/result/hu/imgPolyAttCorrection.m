function imgAttOut =  imgPolyAttCorrection( imgAttIn, spectrum, map )
%   Correct image to standard intensity (Polychromatic version)
%   input:
%       imgAtt 
%       spectrum
%       map
%   output:
%       imgHu
%
% Meng Wu @ stanford
% 2012


[E_soft, mu_soft] = readMaterialAttenuation('Tissue_Soft_ICRU-44', 'physicsdata/materials/');

muSoft = spectralAttenuation(spectrum.energyBinLabels, spectrum.photonsPerEnergyBin, E_soft, mu_soft);

imgTissue = imgAttIn( map.mapTissue ); 

meanTissue = mean( imgTissue(:) );


imgAttOut = imgAttIn * muSoft / meanTissue;



end