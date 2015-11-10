function [E, mu, rho, mu_by_rho, mu_by_rho_en] = readMaterialAttenuation(materialName, materialsDir)
% [E, mu, rho, mu_by_rho, mu_by_rho_en] = readMaterialAttenuation(materialName, materialsDir)
% reads energy-dependent attenuation values from ASCII files of data that
% has been converted from the NIST tables at http://www.nist.gov/pml/data/xraycoef/
%
% INPUT:
% materialName  a string with the material's name (e.g. 'Gold')
% materialsDir  directory with files 'density_*.txt' (containing a single
%   ASCII number with the material's density in g/cm^3) and
%   attenuation_*.txt' (containing the energy levels in keV, the mass
%   attenuation coefficient mu_by_rho in cm^2/g, and the
%   mass energy-absorption coefficient my_by_rho_en in cm^2/g)
%
% OUTPUT:
% E  the energy levels (in keV) at which attenuation values are given
% mu  the attenuation coefficients (in 1/cm) obtained using the standard
%   density rho
% rho  the material's standard density (in g/cm^3)
% mu_by_rho  the corresponding mass attenuation coefficients (in cm^2/g)
% mu_by_rho_en  the corresponding mass energy-absorption coefficients
%   (in cm^2/g)
%
% Copyright (c) 2011 by Andreas Keil, Stanford University.


%% Input Argument Check

if ~ischar(materialName) || ~isvector(materialName)
	error('"materialName" has to be a material''s name!');
end
if ~ischar(materialsDir) || ~isvector(materialsDir) || ~exist(materialsDir, 'dir')
	error('"materialsDir" has to be a string representing an existing directory!');
end


%% Read Densities and Attenuation Tables from Disk

% read density from file
fid = fopen(fullfile(materialsDir, ['density_' materialName '.txt']), 'r');
rho = fscanf(fid, '%f');
fclose(fid);

% read attenuation data from file
fid = fopen(fullfile(materialsDir, ['attenuation_' materialName '.txt']), 'r');
data = fscanf(fid, '%f %f %f', [3 inf]);
fclose(fid);


%% Store Return Values

E = data(1, :);
% make energy values unique around absorption edges
for i = 1:length(E)-1
	if E(i) == E(i+1)
		E(i) = E(i)/(1+2*eps);
		E(i+1) = E(i+1)*(1+2*eps);
	end
end
mu_by_rho = data(2, :);
mu = mu_by_rho * rho;
mu_by_rho_en = data(3, :);
