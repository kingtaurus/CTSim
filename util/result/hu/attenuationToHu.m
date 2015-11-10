function imgHu = attenuationToHu(imgAtt, spectrumEnergies, spectrumPortions, materialsDir)
% imgHu = attenuationToHu(imgAtt, spectrumEnergies, spectrumPortions, materialsDir)
%   computes the Hounsfield units for a given attenuation image; works with
%   any number of image dimensions
%
% INPUT:
% imgAtt  n-D image (or scalar) containing attenuation values
% spectrumEnergies  energies bin of the incident X-ray spectrum (in keV)
% spectrumPortions  portions of the bins used for computing a weighted
%   average attenuation (in arbitraty units)
% materialsDir  path string to the directory containing the material data files
%
% Copyright (c) 2010 by Andreas Keil, Stanford University.


%% Read Attenuation Coefficients of Air and Water

[E_air, mu_air] = readMaterialAttenuation('Air_Dry_near_sea_level', materialsDir);
[E_water, mu_water] = readMaterialAttenuation('Water_Liquid', materialsDir);


%% Generate HU Volume

% compute attenuation of air and water at current energy; material 1 = air, material 2 = water
muAir = spectralAttenuation(spectrumEnergies, spectrumPortions, E_air, mu_air);
muWater = spectralAttenuation(spectrumEnergies, spectrumPortions, E_water, mu_water);

% convert to HU
imgHu = (imgAtt - muWater) / (muWater - muAir) * 1000;
