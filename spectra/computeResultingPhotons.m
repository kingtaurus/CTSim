function photonsOut = computeResultingPhotons(photonsIn, energies, varargin)
% photonsOut = computeResultingPhotons(photonsIn, energies, materialsDir, M1, T1, M2, T2, ...)
%    computes the number of pohtons remaining after traversing a given set
%    of objects for the given vector of photons per energy bin.
%
% INPUT:
% photonsIn   input number of photons, can be a vector
% energies   vector of energies at which the total attenuation is to be
%    computed (keV)
% Mi, Ti   pairs of material names and thicknesses (mm), through which the
%    total attenuation is to be computed
%
% OUTPUT:
% photonsOut   output number of photons after passing through the given
%    materials
% 
% Copyright (c) 2012 by Andreas Keil, Stanford University.


%% Parse Arguments

validateattributes(photonsIn, {'numeric'}, { 'real', 'nonnegative'});
noMaterials = (nargin-2)/2;
materialNames = varargin(1:2:end);
materialThicknesses = varargin(2:2:end); % in mm

%% Compute Resulting Number of Photons
totalAttenuations = zeros(size(energies));
for m = 1:noMaterials
    
	mu = materialAttenuation( energies, materialNames{m});
	totalAttenuations(:) = totalAttenuations(:) + mu(:) * materialThicknesses{m}/10; % thicknesses are given in mm but attenuation values in cm

end

photonsOut = photonsIn .* exp(-totalAttenuations);
