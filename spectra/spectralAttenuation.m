function mu = spectralAttenuation(spectrumEnergies, spectrumPortions, attenuationEnergies, attenuationValues)
% mu = spectralAttenuation(spectrumEnergies, spectrumPortions, attenuationEnergies, attenuationValues)
%   computes the spectrum-specific attenuation using given energy-dependent
%   attenuation values.
%
% INPUT:
% spectrumEnergies  vector of the incident spectrum's energy levels / bin
%   labels (in keV)
% spectrumPortions  vector of the number of photons or ratio per energy bin
%   given in spectrumEnergies
% attenuationEnergies  vector of energy levels at which the input material's
%   attenuation values are given (in keV)
% attenuationValues  vector of the input material's attenuation values at
%   the energy levels given in attenuationEnergies (in 1/cm)
%
% OUTPUT:
% mu  the attenuation coefficients (in 1/cm)
%
% Copyright (c) 2011 by Andreas Keil, Stanford University.


%% Input Argument Check

%TODO: remove checks?
if ~isvector(spectrumEnergies) || ~isvector(spectrumPortions) || length(spectrumEnergies) ~= length(spectrumPortions) || length(spectrumEnergies) < 1
	error('"spectrumEnergies" and "spectrumPortions" have to be non-empty vectors of the same length!');
end
if ~isvector(attenuationEnergies) || ~isvector(attenuationValues) || length(attenuationEnergies) ~= length(attenuationValues) || length(attenuationEnergies) < 1
	error('"attenuationEnergies" and "spectrumPortions" have to be non-empty vectors of the same length!');
end


%% Compute Spectrum-Specific Attenuation as Weighted Average of Energy-Dependent Attenuation Values

% mu = 0;
% for i = 1:length(spectrumEnergies)
% 	mu = mu + spectrumPortions(i) * interp1geom(attenuationEnergies, attenuationValues, spectrumEnergies(i));
% end
% mu = mu / sum(spectrumPortions);
mu = spectralAttenuations(spectrumEnergies, attenuationEnergies, attenuationValues);
mu = dot(spectrumPortions, mu) / sum(spectrumPortions);
