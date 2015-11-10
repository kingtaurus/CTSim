function p = readParametersKVMV(filename)
% p = readParameters(filename) 
%   reads the parameters defined in this function from the given ini file
%   and returns them in a struct
%
% Copyright (c) 2012 by Andreas Keil, Stanford University.


% define parameters
parameterDefs = { ...
	'Phantom',        '',       'materialsFileName',                     's'; ... % phantom file with material indexes
	'Phantom',        '',       'materialMappingName',                   's'; ... % name of the material mapping to apply
	...
	'Reconstruction', '',       'size',                                  'd'; ... % reconstruction size (vector of pixels)
	'Reconstruction', '',       'spacing',                               'd'; ... % reconstruction spacing (vector of mm)
    'Reconstruction', '',       'offset',                                'd'; ... % reconstruction spacing (vector of mm)
	...
	'Spectra',        '',       'spectrumKeV',                           's'; ... % keV spectrum to use
	'Spectra',        '',       'spectrumMeV',                           's'; ... % MeV spectrum to use
	'Spectra',        '',       'doseLimitTotalKeV',                     'd'; ... % total dose limit for keV imaging (Gy)
	'Spectra',        '',       'doseLimitTotalMeV',                     'd'; ... % total dose limit for MeV imaging (Gy)
	...
	'Geometries',     'keV',    'SAD',                                   'd'; ... % source-to-axis distance (mm)
	'Geometries',     'keV',    'ADD',                                   'd'; ... % axis-to-detector distance (mm)
	'Geometries',     'keV',    'noViews',                               'd'; ... % number of views/projections
	'Geometries',     'keV',    'sizeDet',                               'd'; ... % number of detector pixels (per acqusition angle / views)
	'Geometries',     'keV',    'spacingDet',                            'd'; ... % detector pixel spacing (mm)
    'Geometries',     'keV',    'offsetDet',                             'd'; ... % detector pixel spacing (mm)
    'Geometries',     'keV',    'sinoBasedSegm',                         'd'; ... % metal segmentation method
	'Geometries',     'MeV',    'SAD',                                   'd'; ... % source-to-axis distance (mm)
	'Geometries',     'MeV',    'ADD',                                   'd'; ... % axis-to-detector distance (mm)
	'Geometries',     'MeV',    'noViews',                               'd'; ... % number of views/projections
	'Geometries',     'MeV',    'sizeDet',                               'd'; ... % number of detector pixels (per acqusition angle / views)
	'Geometries',     'MeV',    'spacingDet',                            'd'; ... % detector pixel spacing (mm)
    'Geometries',     'MeV',    'offsetDet',                             'd'; ... % detector pixel spacing (mm)
	...
	'Detector',       '',       'detectorConversionEfficiencyKeV',       'd'; ... % fraction of photons that are actually detected/converted by the detector
	'Detector',       '',       'detectorConversionEfficiencyMeV',       'd'; ... % fraction of photons that are actually detected/converted by the detector
	...
	'Paths',          '',       'materialsDir',                          's'; ... % path to the directory with the material data (mass attenuation coefficients and density)
	'Paths',          '',       'spectraDir',                            's'; ... % path to the directory with the spectra tables
	...
	'Visualization',  '',       'windowAtt',                             'd'; ... % standard window for displaying attenuation values
	'Visualization',  '',       'windowHu',                              'd'; ... % standard window for displaying Hounsfield units
	'Visualization',  '',       'windowSinoKeV',                         'd'; ... % standard window for displaying attenuation integral sinograms with keV energies
	'Visualization',  '',       'windowSinoMeV',                         'd'; ... % standard window for displaying attenuation integral sinograms with keV energies
	...
	'Bowtie',           '',       'shapeType',                           's'; 
	'Bowtie',           '',       'alpha',                               'd'; 
	'Bowtie',           '',       'beta',                                'd'; 
	'Bowtie',           '',       'maximumThickness',                    'd'; 
    'Bowtie',           '',       'minimumThickness',                    'd'; 
    'Bowtie',           '',       'material',                            's'; 
    };
noParameters = length(parameterDefs);

% read parameters from file
[values, errorMsgs] = inifile(filename, 'read', parameterDefs);
if length(values) ~= noParameters || length(errorMsgs) ~= noParameters
	error('Error reading parameters!  Returned cell arrays should have same length as parameter definitions.');
end
anyErrors = false;
for i = 1:length(errorMsgs)
	if ~isempty(errorMsgs{i})
		anyErrors = true;
		disp(errorMsgs{i});
	end
end
if anyErrors
	error('Error reading parameters!');
end

% populate workspace with variables
for i = 1:noParameters
	variableGroup = strrep(parameterDefs{i, 1}, ' ', '');
	variableSubGroup = strrep(parameterDefs{i, 2}, ' ', '');
	variableName = strrep(parameterDefs{i, 3}, ' ', '');
	fullVariableName = 'p';
	if ~isempty(variableGroup)
		fullVariableName = [fullVariableName '.' variableGroup]; %#ok<AGROW>
	end
	if ~isempty(variableSubGroup)
		fullVariableName = [fullVariableName '.' variableSubGroup]; %#ok<AGROW>
	end
	fullVariableName = [fullVariableName '.' variableName]; %#ok<AGROW>
	eval([fullVariableName '= values{i};']);
end




end


