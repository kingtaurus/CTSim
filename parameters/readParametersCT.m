function p = readParametersCT(filename)
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
    'Spectra',        '',       'spectrum',                           's'; ... % keV spectrum to use
    'Spectra',        '',       'focalSpotSize',                      'd'; ... % focal spot size of the x-ray tube
    'Spectra',        '',       'maximumIntensity',                   'd'; ... % number of photon counts in mm^2 at 1000 mm away
    'Spectra',        '',       'automaticExposureControl',           'd'; ... % enable mA modulation for AEC
    ...
    'Geometries',     '',    'SAD',                                   'd'; ... % source-to-axis distance (mm)
    'Geometries',     '',    'ADD',                                   'd'; ... % axis-to-detector distance (mm)
    'Geometries',     '',    'noViews',                               'd'; ... % number of views/projections
    'Geometries',     '',    'sizeDet',                               'd'; ... % number of detector pixels (per acqusition angle / views)
    'Geometries',     '',    'spacingDet',                            'd'; ... % detector pixel spacing (mm)
    'Geometries',     '',    'offsetDet',                             'd'; ... % detector pixel spacing (mm)
    'Geometries',     '',    'flatPanel',                             'd'; ...
    ...
    'Detector',       '',       'detectorConversionEfficiency',         'd'; ... % fraction of photons that are actually detected/converted by the detector
    'Detector',       '',       'pointSpreadFunctionFWHM',              'd'; ... % fraction of photons that are actually detected/converted by the detector
    'Detector',       '',       'noisePowerSpectrum',                   'd'; ... % fraction of photons that are actually detected/converted by the detector
    'Detector',       '',       'energyIntegrating',                    'd'; ...
    'Detector',       '',       'compoundPoissonNoise',                 'd'; ...
    ...
    'Paths',          '',       'materialsDir',                          's'; ... % path to the directory with the material data (mass attenuation coefficients and density)
    'Paths',          '',       'spectraDir',                            's'; ... % path to the directory with the spectra tables
    ...
    'Visualization',  '',       'windowAtt',                             'd'; ... % standard window for displaying attenuation values
    'Visualization',  '',       'windowHu',                              'd'; ... % standard window for displaying Hounsfield units
    'Visualization',  '',       'windowSinoKeV',                         'd'; ... % standard window for displaying attenuation integral sinograms with keV energies
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


