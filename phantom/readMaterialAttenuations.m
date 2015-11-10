function [E, mu, rho, mu_by_rho, mu_by_rho_en] = readMaterialAttenuations(materials, materialsDir)
% [E, mu, rho, mu_by_rho, mu_by_rho_en] = readMaterialAttenuations(materials, materialsDir)
%   reads the densities and attenuation tables for the given materials
%
% INPUT:
% materials  a cell array of matrial names (the corresponding attenuation
%   data is read from files in the matrialsDir)
% materialsDir  path string to the directory containing the material data files
%
% OUTPUT:
% E  cell (one per material) of vectors of photon energy values (in keV) for which
%   for which the attenuation values are returned in mu_by_rho and mu
% mu  cell (one per material) of vectors of attenuation values (in 1/cm)
%   corresponding to the energies in E and the standard density rho
% rho  cell of densities of the given materials
% mu_by_rho  cell (one per material) of vectors of mass attenuation
%   coefficients (in cm^2/g) corresponding to the energies in E
% mu_by_rho_en  cell (one per material) of vectors of mass energy-absorption
%   cofficients (in cm^2/g) corresponding to the energies in E
%
% Copyright (c) 2011 by Andreas Keil, Stanford University.


%% Input Argument Check

if ~iscell(materials) || ~isvector(materials)
	error('"materials" has to be a cell vector!');
end


%% Read Material-Dependent Attenuation Coefficients

E = cell(size(materials));
mu = cell(size(materials));
rho = cell(size(materials));
mu_by_rho = cell(size(materials));
mu_by_rho_en = cell(size(materials));
for m = 1:length(materials)
	[E{m}, mu{m}, rho{m}, mu_by_rho{m}, mu_by_rho_en{m}] = readMaterialAttenuation(materials{m}, materialsDir);
end
