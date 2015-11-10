function mu = spectralAttenuations(spectrumEnergies, attenuationEnergies, attenuationValues)
% mu = spectralAttenuations(spectrumEnergies, attenuationEnergies, attenuationValues)
%   computes the spectrum-specific attenuation using given energy-dependent
%   attenuation values.
%
% INPUT:
% spectrumEnergies  vector of the incident spectrum's energy levels / bin
%   labels (in keV)
% attenuationEnergies  vector of energy levels at which the input material's
%   attenuation values are given (in keV)
% attenuationValues  vector of the input material's attenuation values at
%   the energy levels given in attenuationEnergies (in 1/cm)
%
% OUTPUT:
% mu  the attenuation coefficients at the given spectrumEnergies (in 1/cm)
%
% Copyright (c) 2012 by Andreas Keil, Stanford University.


%% Input Argument Check

%TODO: remove checks?
if ~isvector(spectrumEnergies) || length(spectrumEnergies) < 1
	error('"spectrumEnergies" has to be a non-empty vector!');
end
if ~isvector(attenuationEnergies) || ~isvector(attenuationValues) || length(attenuationEnergies) ~= length(attenuationValues) || length(attenuationEnergies) < 1
	error('"attenuationEnergies" and "attenuationValues" have to be non-empty vectors of the same length!');
end


%% Compute Spectrum-Specific Attenuation as Weighted Average of Energy-Dependent Attenuation Values

mu = zeros(length(spectrumEnergies), 1);
for i = 1:length(spectrumEnergies)
	mu(i) = interp1geom(attenuationEnergies, attenuationValues, spectrumEnergies(i));
end
