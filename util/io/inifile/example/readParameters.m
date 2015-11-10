% readParameters is a script that reads the parameters defined below from
% 'parameters.ini' in the current path and stores them in the workspace
%
% Copyright (c) 2012 by Andreas Keil, Stanford University.


% define parameters
parameterDefs = { ...
	'Spectra',                        '', 'spectrumKeV',                     's'; ... % keV spectrum to use
	'Spectra',                        '', 'spectrumMeV',                     's'; ... % MeV spectrum to use
	...
	'Algorithm Parameters',           '', 'impactCoefficientsFileName',      's'; ... % the mu(70keV) -> Phi/Theta mapping to use for Bruno De Man's IMPACT algorithm
	'Algorithm Parameters',           '', 'impactE0',                        'd'; ... % reference energy
	'Algorithm Parameters',           '', 'beamHardeningCorrect',            'd'; ... % 
	'Algorithm Parameters',           '', 'dummyLowPhotonCount',             'd'  ... % the dummy value that replaces any photon count < 1
	};
noParameters = length(parameterDefs);

% read parameters from file
[values, errorMsgs] = inifile('parameters.ini', 'read', parameterDefs);
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
	variableName = parameterDefs{i, 3};
	eval([variableName, '= values{i};']);
end

% clean up everything that was not a parameter
clear parameterDefs noParameters values errorMsgs anyErrors variableName i
