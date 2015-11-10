function [energyBinLabels, energyBinWidths, photonsPerEnergyBin, photonDensity, photonsTotal, energyAverage] = readSpectrum(filename, densityOrPerBin, noPhotonsScaleFactor)
% [energyBinLabels,
%  energyBinWidths,
%  photonsPerEnergyBin,
%  photonDensity,
%  photonsTotal,
%  energyAverage]
% = readSpectrum(filename, densityOrPerBin, noPhotonsScaleFactor)
%
% reads an X-ray spectrum from file and returns the data in various formats.
%
% INPUT:
% filename   name of a file containing two columns of values in ASCII
%   format separated by a space, where column 1 contains the energy labels
%   and column 2 the number of photons per mm^2 in 1m distance from the
%   source
% densityOrPerBin   if 'per bin' then the number of photons are given per
%   bin (where the boundaries of each bin are half-way between the bin
%   labels);
%   if 'relative' then all spectrum data is normalized to sum up to 1 and
%   then multiplied by noPhotonsScaleFactor (which is the total number of
%   photons in this case) to yield the number of photons per bin;
%   if 'density' then the number of photons are samples of the density
%   function per keV
% noPhotonsScaleFactor   a scalar factor (e.g., depending on the tube's mA
%   setting, acquisition time, source-to-detector distance, and pixel size)
%   that is used to scale all number of photons
%
% OUTPUT:
% energyBinLabels   copies of the energy levels given in column 1 of the
%   data file
% energyBinWidths   the widths of the energy bins (where the boundaries of
%   each bin are half-way between the bin labels)
% photonsPerEnergyBin   the number of photons per bin
% photonDensity   the sampled photon density function (per keV)
% photonsTotal   the total number of photons of all energies
% energyAverage   the weighted average energy of the photons
%
% Note all data rows are stripped of the zero photon rows.
%
% Copyright (c) 2011 by Andreas Keil, Stanford University.


%% Input Check

if ~exist(filename, 'file')
	error('Spectrum file "%s" does not exist!', filename);
end
if ~isscalar(noPhotonsScaleFactor) || ~isreal(noPhotonsScaleFactor) || noPhotonsScaleFactor <= 0
	error('noPhotonsScaleFactor has to be a real, positive scalar!');
end


%% Debug Output

%fprintf(['Reading spectrum from "' filename '"... ']);


%% Read File

fid = fopen(filename, 'r');
spectrum = fscanf(fid, '%f %f', [2 inf]); spectrum = spectrum';
fclose(fid);


%% Pre-Process Raw Data

% check that at least first and last bin have no photons
if (spectrum(1, 2) ~= 0) || (spectrum(end, 2) ~= 0)
	error('First and last energy bin in the given file must be empty (with 0 photons)!');
end

% remove all data rows with no photons (except the first and the last one)
for firstNonZero = 2:size(spectrum, 1)-1
	if spectrum(firstNonZero, 2) ~= 0
		break;
	end
end
for lastNonZero = size(spectrum, 1)-1:-1:2
	if spectrum(lastNonZero, 2) ~= 0
		break;
	end
end
spectrum = spectrum(firstNonZero-1:lastNonZero+1, :);


%% Store/Compute Output Variables

energyBinLabels = spectrum(:, 1);
if strcmpi(densityOrPerBin, 'per bin')
	% scale the number of photons by the given factor
	photonsPerEnergyBin = spectrum(:, 2) * noPhotonsScaleFactor;
	% compute # of photons per keV using the centers between bins as bin boundaries
	energyBinWidths = (energyBinLabels(3:end)-energyBinLabels(1:end-2)) / 2;
	photonDensity = photonsPerEnergyBin(2:end-1) ./ energyBinWidths;
	% discard the first and last bin since they contain no photons
	energyBinLabels = energyBinLabels(2:end-1);
	photonsPerEnergyBin = photonsPerEnergyBin(2:end-1);
elseif strcmpi(densityOrPerBin, 'relative')
	warning('Double-check this portion of the code.  It''s not been tested thouroughly.');
	% compute relative # of photons per bin and then multiply with total number of photons given
	photonsPerEnergyBin = spectrum(:, 2)/sum(spectrum(:, 2)) * noPhotonsScaleFactor;
	% compute # of photons per keV using the centers between bins as bin boundaries
	energyBinWidths = (energyBinLabels(3:end)-energyBinLabels(1:end-2)) / 2;
	photonDensity = photonsPerEnergyBin(2:end-1) ./ energyBinWidths;
	% discard the first and last bin since they contain no photons
	energyBinLabels = energyBinLabels(2:end-1);
	photonsPerEnergyBin = photonsPerEnergyBin(2:end-1);
elseif strcmpi(densityOrPerBin, 'density')
	warning('Double-check this portion of the code.  It''s not been tested thouroughly.');
	% store # of photons per keV, discarding the first and last bin since they contain no photons
	photonDensity = spectrum(2:end-1, 2);
	% compute # of photons per bin using the centers between bins as bin boundaries
	energyBinWidths = (energyBinLabels(3:end)-energyBinLabels(1:end-2)) / 2;
	photonsPerEnergyBin = photonDensity .* energyBinWidths * noPhotonsScaleFactor;
	% discard the first and last bin since they contain no photons
	energyBinLabels = energyBinLabels(2:end-1);
else
	error('Unknown input value ''%s'' for densityOrPerBin argument! Use either ''denisty'' or ''per bin''.');
end

% compute the total number of photons
photonsTotal = sum(photonsPerEnergyBin);

% compute average photon energy
energyAverage = energyBinLabels' * photonsPerEnergyBin / photonsTotal;

% some quick final tests
photonsTotalFromDensity = photonDensity' * energyBinWidths;
if abs(photonsTotal - photonsTotalFromDensity)/photonsTotal > 100*eps
	error('%f ~= %f !', photonsTotal, photonsTotalFromDensity);
end
photonsTotalFromPerBinValues = sum(photonsPerEnergyBin);
if abs(photonsTotal - photonsTotalFromPerBinValues)/photonsTotal > 100*eps
	error('%f ~= %f !', photonsTotal, photonsTotalFromPerBinValues);
end


%% Debug Output

%fprintf('done.\n');

