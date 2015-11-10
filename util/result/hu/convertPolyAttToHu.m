function imgHu = convertPolyAttToHu( imgAtt, spectrum)
%   Convert attenuation to Hounsfield Units (Polychromatic version)
%   input:
%       imgAtt 
%       spectrum
%   output:
%       imgHu
%
% Meng Wu @ stanford
% 2012


% Generate HU Volume
[~, muWater] = materialAttenuation(spectrum.energyBinLabels, 'water', spectrum.photonsPerEnergyBin);

% convert to HU
imgHu = (imgAtt - muWater) / (muWater) * 1000;


end