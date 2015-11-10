function imgAttPoly = impactImageToPolychromaticAttenuation(imgAttE0, E0, deManMu, deManPhi, deManTheta, deManKappa, deManPi, modelEdge, EKappa, epsilonKappa, modelPairProduction, spectrumEnergies, spectrumPortions)
% imgAtt = impactImageToPolychromaticAttenuation(imgAttE0, ...
%    E0, deManMu, deManPhi, deManTheta, deManKappa, deManPi, ...
%    modelEdge, EKappa, epsilonKappa, modelPairProduction, ...
%    spectrumEnergies, spectrumPortions)
% returns the attenuation image for a given image of material indexes and
% at a given photon energy; works with any number of image dimensions
%
% INPUT:
% imgAttE0
% E0, deManMu, deManPhi, deManTheta, deManKappa, deManPi
% modelEdge, EKappa, epsilonKappa, modelPairProduction
% spectrumEnergies   vector of the incident spectrum's energy levels / bin labels
% spectrumPortions   vector of the number of photons or ratio per energy bin
%   given in spectrumEnergies
%
% OUTPUT:
% imgAttPoly   n-D image of polychromatic attenuation values
%
% Copyright (c) 2012 by Andreas Keil, Stanford University.


imgPhi = interp1(deManMu, deManPhi, imgAttE0, 'linear', 'extrap');
imgTheta = interp1(deManMu, deManTheta, imgAttE0, 'linear', 'extrap');
if modelEdge
	imgKappa = interp1(deManMu, deManKappa, imgAttE0, 'linear', 'extrap');
end
if modelPairProduction
	imgPi = interp1(deManMu, deManPi, imgAttE0, 'linear', 'extrap');
end

imgAttPoly = zeros(size(imgAttE0));
for i = 1:length(spectrumEnergies)
	Ei = spectrumEnergies(i);
	imgAttEi = imgPhi * photoElectric(Ei)/photoElectric(E0) + imgTheta * kleinNishina(Ei)/kleinNishina(E0);
	if modelEdge, imgAttEi = imgAttEi + imgKappa * photoElectric(Ei)/photoElectric(E0) * HKappa(Ei, EKappa, epsilonKappa); end
	if modelPairProduction, imgAttEi = imgAttEi + imgPi * pairProduction(Ei); end
	imgAttPoly = imgAttPoly + spectrumPortions(i)*imgAttEi;
end
imgAttPoly = imgAttPoly / sum(spectrumPortions);
