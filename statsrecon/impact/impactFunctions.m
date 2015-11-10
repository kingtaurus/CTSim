function [muE0, Phi, Theta, Kappa, Pi] = impactFunctions(...
	materialNames, ...
	materialsDir, ...
	E0, ...
	Emin, ...
	Emax, ...
	modelEdge, ...
	EKappa, ...
	epsilonKappa, ...
	modelPairProduction, ...
	interactive)
% [muE0, Phi, Theta, Kappa, Pi] = impactFunctions(...
%    materialNames, ...
%    materialsDir, ...
%    E0, ...
%    Emin, ...
%    Emax, ...
%    modelEdge, ...
%    EKappa, ...
%    epsilonKappa, ...
%    modelPairProduction, ...
%    interactive)
%
%    generates the parameter lists of coefficients Phi, Theta, Kappa, and Pi
%    from Bruno De Man's IMPACT Algorithm (Trans. Med. Imag. 20(10), pp.
%    999-1008, 2001) modelling the material dependence of X-ray attenuation
%    for a set of given materials.
%
% INPUT:
% materialNames   
% materialsDir   
% E0   
% Emin   
% Emax   
% modelEdge   should a third coefficient be added for K-edge modelling?
%    (true/false)
% EKappa   
% epsilonKappa   
% modelPairProduction   should a forth coefficient be added for modelling
%    pair production? (true/false)
% interactive   should fitted attenuation profiles be plotted on screen?
%    (true/false, default = true)
%
% OUTPUT:
% muE0    attenuation at E0
% Phi     De Man coefficient for photoelectric absorption
% Theta   De Man coefficient for Compton scattering
% Kappa   Andreas Keil's additional De Man coefficient for photoelectric
%    absorption above a K edge
% Pi      Andreas Keil's additional De Man coefficient for absorption
%    through pair production above 1.022MeV
%
% Copyright (c) 2012 by Andreas Keil, Stanford University.


%% Parse Input Arguments

validateattributes(materialNames, {'cell'}, {});
validateattributes(materialsDir, {'char'}, {});
validateattributes(E0, {'numeric'}, {'scalar', 'real', 'positive'});
validateattributes(Emin, {'numeric'}, {'scalar', 'real', 'positive'});
validateattributes(Emax, {'numeric'}, {'scalar', 'real', 'positive'});
validateattributes(modelEdge, {'logical'}, {'scalar'});
validateattributes(EKappa, {'numeric'}, {'scalar', 'real', 'positive'});
validateattributes(epsilonKappa, {'numeric'}, {'scalar', 'real', 'positive'});
validateattributes(modelPairProduction, {'logical'}, {'scalar'});
if nargin < 10, interactive = true; end
validateattributes(interactive, {'logical'}, {'scalar'});


%% Check Output Arguments

if nargout < 3 || (modelEdge && nargout < 4) || (modelPairProduction && nargout < 5)
	error('Too few output arguments requested!');
end


%% Define Materials Included in the Overall De Man Model

noMaterials = length(materialNames);
if noMaterials == 1 && strcmp(materialNames{1}, 'all')
	materialNames = dir([materialsDir 'attenuation_*.txt']);
	materialNames = {materialNames.name};
	noMaterials = length(materialNames);
	for i = 1:noMaterials
		materialNames{i} = strrep(materialNames{i}, 'attenuation_', '');
		materialNames{i} = strrep(materialNames{i}, '.txt', '');
	end
end


%% Compute Parameters for Each Material

% initialize parameters
muE0 = zeros(noMaterials, 1);
Phi = zeros(noMaterials, 1);
Theta = zeros(noMaterials, 1);
Kappa = zeros(noMaterials, 1);
Pi = zeros(noMaterials, 1);
edgeEnergy = zeros(noMaterials, 1);

% compute scalar parameters Phi and Theta for each material
for i = 1:noMaterials
	[muE0(i) Phi(i) Theta(i) Kappa(i) Pi(i) edgeEnergy(i)] = ...
		impactFitCoefficients(...
		materialNames{i}, ...
		materialsDir, ...
		E0, ...
		Emin, ...
		Emax, ...
		modelEdge, ...
		EKappa, ...
		epsilonKappa, ...
		modelPairProduction, ...
		interactive);
end

% sort data by attenuation value at E0
[muE0 sortIndex] = sort(muE0);
materialNames = materialNames(sortIndex);
Phi = Phi(sortIndex);
Theta = Theta(sortIndex);
Kappa = Kappa(sortIndex);
Pi = Pi(sortIndex);
edgeEnergy = edgeEnergy(sortIndex);


%% Output List of Coefficients on Screen

fprintf('Material Name                  |mu @%4gkeV|   Phi   |  Theta  |  Kappa  |   Pi    |main edge\n', E0);
fprintf('-------------------------------+-----------+---------+---------+---------+---------+---------\n');
for i = 1:noMaterials
	if edgeEnergy(i) == 0
		edgeStr = '(no edge)';
	else
		edgeStr = sprintf('%.2fkeV', edgeEnergy(i));
	end
	fprintf('%-31s| %10.6f|%9.6f|%9.6f|%9.6f|%9.6f|%s\n', materialNames{i}, muE0(i), Phi(i), Theta(i), Kappa(i), Pi(i), edgeStr);
end


%% Plot Coefficient Functions

if interactive
	% prepare legend labels
	labels = {'Phi', 'Theta'};
	if modelEdge, labels = [labels, {'Kappa'}]; end
	if modelPairProduction, labels = [labels, {'Pi'}]; end
	
	% plot actual attenuation and approximation on linear axes
	figure('Name', sprintf('Mapping Attenuation @ %g keV -> Attenuation Components', E0), 'NumberTitle', 'off');
	zoom on;
	plot(muE0, Phi, 'x-');
	hold all;
	plot(muE0, Theta, 'o-');
	if modelEdge, plot(muE0, Kappa, 's-'); end
	if modelPairProduction, plot(muE0, Pi, '+-'); end
	for i = 1:noMaterials
		text(muE0(i), mean([Phi(i), Theta(i)]), materialNames{i}, 'Interpreter', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Rotation', 90);
	end
	grid on;
	xlabel(sprintf('\\mu_{%gkeV} (1/cm)', E0));
	ylabel('Attenuation component parameters (1/cm)');
	legend(labels, 'Location', 'NorthWest');

	
	% plot actual attenuation and approximation on log log axes
	figure('Name', sprintf('Mapping Attenuation @ %g keV -> Attenuation Components', E0), 'NumberTitle', 'off');
	zoom on;
	loglog(muE0, Phi, 'x-');
	hold all;
	loglog(muE0, Theta, 'o-');
	if modelEdge, loglog(muE0, Kappa, 's-'); end
	if modelPairProduction, loglog(muE0, Pi, '+-'); end
	for i = 1:noMaterials
		text(muE0(i), sqrt(Phi(i)*Theta(i)), materialNames{i}, 'Interpreter', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Rotation', 90);
	end
	grid on;
	xlabel(sprintf('\\mu_{%gkeV} (1/cm)', E0));
	ylabel('Attenuation component parameters (1/cm)');
	legend(labels, 'Location', 'NorthWest');
end
