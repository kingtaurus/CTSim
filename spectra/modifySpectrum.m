function modifySpectrum(materialsDir, filterMaterial, filterThickness, downsamplingFactor)
% modifySpectrum()  is used to filter and downsample a spectrum.  The input
% file name is queried interactively through a dialog and the output file
% name is constructed automatically by appending the filtration and
% downsampling parameters.
%
% INPUT:
% materialsDir   directory where to find the materials data used for
%    filtering spectra
% filterMaterial   name of the filtration material
% filterThickness  thickness of the filter (mm)
% downsamplingFactor   factor by which a spectrum is to be downsampled
%
% Copyright (c) 2012 by Andreas Keil, Stanford University.


%% Check Input Arguments

validateattributes(materialsDir, {'char'}, {});
validateattributes(filterMaterial, {'char'}, {});
if isempty(materialsDir) || isempty(filterMaterial)
	filterSpectrum = false;
else
	validateattributes(filterThickness, {'numeric'}, {'scalar', 'real', 'nonnegative'});
	if filterThickness == 0
		filterSpectrum = false;
	else
		filterSpectrum = true;
	end
end

if nargin < 4
	downsamplingFactor = 1;
end
validateattributes(downsamplingFactor, {'numeric'}, {'scalar', 'integer', 'positive'});
downsampleSpectrum = (downsamplingFactor ~= 1);


%% Query File Names

% ask for input file name
filters = { ...
	'*.txt',  'Text files (*.txt)'; ...
	'*.*',  'All Files (*.*)' ...
	};
[file, path] = uigetfile(filters, 'Load Original Spectrum');
if isscalar(file) && isscalar(path) && isequal(file, 0) && isequal(path, 0)
	disp('No filename specified - aborting.');
	return;
else
	filenameIn = [path, file];
end

% construct output file name
[path, name, ext] = fileparts(filenameIn);
if filterSpectrum
	filterStr = sprintf('.filter_%gmm_%s', filterThickness, filterMaterial);
else
	filterStr = '';
end
if downsampleSpectrum
	downsampleStr = sprintf('.down_%i', downsamplingFactor);
else
	downsampleStr = '';
end
filenameOut = [path, '/', name, filterStr, downsampleStr, ext];


%% Convert Spectrum

% read input spectrum
fidIn = fopen(filenameIn, 'r');
spectrumIn = fscanf(fidIn, '%f %f', [2 inf]); spectrumIn = spectrumIn';
fclose(fidIn);

% just copy top and bottom lines with zero photons
spectrumOutTop = zeros(0, 2);
while spectrumIn(1, 2) == 0
	spectrumOutTop = [spectrumOutTop; spectrumIn(1, :)];
	spectrumIn = spectrumIn(2:end, :);
end
spectrumOutBottom = zeros(0, 2);
while spectrumIn(end, 2) == 0
	spectrumOutBottom = [spectrumOutBottom; spectrumIn(end, :)];
	spectrumIn = spectrumIn(1:end-1, :);
end

% filter spectrum
if filterSpectrum
	[E, mu] = readMaterialAttenuation(filterMaterial, materialsDir);
	attenuations = interp1(E, mu, spectrumIn(:, 1));
	attenuations = attenuations*filterThickness/10;
	spectrumIn(:, 2) = spectrumIn(:, 2) .* exp(-attenuations);
end

% downsample spectrum
if downsampleSpectrum
	noEnergiesOut = ceil(size(spectrumIn, 1) / downsamplingFactor);
	spectrumOut = zeros(noEnergiesOut, 2);
	for i = 1:noEnergiesOut
		firstIndex = (i-1)*downsamplingFactor+1;
		lastIndex = min(i*downsamplingFactor, size(spectrumIn, 1));
		energiesIn = spectrumIn(firstIndex:lastIndex, 1);
		photonsIn = spectrumIn(firstIndex:lastIndex, 2);
		photonsOut = sum(photonsIn);
		energyOut = (energiesIn' * photonsIn) / photonsOut;
		spectrumOut(i, 1) = energyOut;
		spectrumOut(i, 2) = photonsOut;
	end
else
	spectrumOut = spectrumIn;
end

% concatenate downsampled spectrum with zero lines
spectrumOut = [spectrumOutTop; spectrumOut; spectrumOutBottom];

% write output spectrum
fidOut = fopen(filenameOut, 'w');
fprintf(fidOut, '%10.3e %14.7e\n', spectrumOut');
fclose(fidOut);
