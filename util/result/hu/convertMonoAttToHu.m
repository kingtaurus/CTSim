function imgHu = convertMonoAttToHu( imgAtt, spectrum)
% Convert attenuation coefficients to Hounsfield units ( Monochromatic version)
% input:   
%       imgAtt - attenuation coefficients
%       spectrum - x-ray spectrum
% output:
%       imgHu - Hounsfield units
%
%
% Based on Andreas Keil's code       
% 2012.9 Meng Wu at Stanford University 

muWater = materialAttenuation(spectrum.energyAverage, 'water');

% convert to HU
imgHu = (imgAtt - muWater) / muWater  * 1000;

end