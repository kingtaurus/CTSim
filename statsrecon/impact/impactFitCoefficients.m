function [muE0, Phi, Theta, Kappa, Pi, edgeEnergy] = ...
	impactFitCoefficients(materialName, materialsDir, ...
	E0, Emin, Emax, ...
	modelEdge, EKappa, epsilonKappa, ...
	modelPairProduction, ...
	interactive)
% [muE0 Phi Theta Kappa Pi edgeEnergy]
% = impactFitCoefficients(materialName, materialsDir,
% E0, Emin, Emax,
% modelEdge, EKappa, epsilonKappa,
% modelPairProduction,
% interactive)
%
% computes the two scalar parameters Phi and Theta that are used to model
% the material dependence of the given material's attenuation profile in
% Bruno De Man's IMPACT Algorithm (Trans. Med. Imag. 20(10), pp. 999-1008, 2001)
%
% INPUT:
% materialName   name of the material to compute the two parameters for
% materialsDir   directory where the NIST tables are stored
% E0   reference energy (keV)
% Emin   lower limit of energies in the material's attenuation profile that
%    are considered in the approximation (keV)
% Emax   lower limit of energies in the material's attenuation profile that
%    are considered in the approximation (keV)
% (neglecting low and high energies in the approximation helps to focus on
% the two dominant effects that are relevant in X-ray imaging)
% modelEdge   boolean value determining whether a K-edge should be
%    included in the model or not
% EKappa   energy of the one K-edge modeled by our modified IMPACT
%    functions (keV)
% epsilonKappa   scalar value determining the smoothness of the modeled
%    edge at EKappa (keV)
% modelPairProduction   boolean value determining whether pair production
%    should be included in the model or not
% interactive   whether or not to plot (default = true)
%
% OUTPUT:
% muE0    attenutation at the reference energy E0 (1/cm)
% Phi     scalar parameter Phi that is used in the photoelectric part of
%    De Man's attenuation model (1/cm)
% Theta   scalar parameter Theta that is used in the Comption scattering
%    part of De Man's attenuation model (1/cm)
%
% OPTIONAL OUTPUT:
% Kappa   scalar parameter Kappa that is used to model the additional
%    photoelectric attenuation above a K-edge (1/cm)
% Pi   scalar parameter Pi that is used to model the attenuation through
%    pair production (1/cm)
% edgeEnergy   energy of the given material's biggest edge that is not
%    covered by De Man's two-component model (keV)
%
% Copyright (c) 2012 by Andreas Keil, Stanford University.



%% Parse Input Arguments

validateattributes(materialName, {'char'}, {});
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


%% Get Material's Attenuation Profile

% read attenuation coefficients
[E, mu] = readMaterialAttenuation(materialName, materialsDir);

% make sure we have column vectors
E = E(:);
mu = mu(:);


%% Add Some Relevant Energy Values

% add attenuations at E0 and around EKappa
EMeV = 1200;
muE0 = interp1geom(E, mu, E0); % compute this before adding other values just to make sure the interpolation works with unique input data points
mu1MeV = interp1geom(E, mu, EMeV); % compute this before adding other values just to make sure the interpolation works with unique input data points
EKappaRange = (EKappa-epsilonKappa:epsilonKappa/10:EKappa+epsilonKappa)';
muEKappaRange = interp1geom(E, mu, EKappaRange);
E = [E; E0; EKappaRange; EMeV];
mu = [mu; muE0; muEKappaRange; mu1MeV];

% cut-off low (< Emin) and high (> Emax) energies
retainData = (E >= Emin) & (E <= Emax);
E = E(retainData);
mu = mu(retainData);

% make sure data values are sorted and unique
[E, idx] = sort(E);
mu = mu(idx);
[E, idx] = unique(E);
mu = mu(idx);


%% Identify Largest Edge, i.e., the Largest Relative Step in Attenuation that Could Lead to a Large Approximation Error

edgeEnergy = 0;
muRatioMax = 0;
for i = 1:length(E)-1
	if abs(E(i+1)/E(i)-1) < 7*eps % edges are represented in the NIST tables by specifying two attenutation values for two energy values that are very close to each other
		if mu(i+1)/mu(i) > muRatioMax % look at the relative change in attenuation for the edge's size and retain the biggest edge's information
			edgeEnergy = mean(E(i:i+1));
			muRatioMax = mu(i+1)/mu(i);
		end
	end
end


%% Optimize the Material Parameters (for Photoelectric and Compton Effects) for this Material

% determine whether or not we should fit a third, new IMPACT coefficient
fitEdge = modelEdge && ...
	edgeEnergy ~= 0 && ...
	Emin < edgeEnergy && edgeEnergy < Emax && ...
	Emin < EKappa+epsilonKappa && EKappa-epsilonKappa < Emax;

% set up the "right side" and the matrix of the linear least squares problem
b = mu;
A = [photoElectric(E)./photoElectric(E0), kleinNishina(E)./kleinNishina(E0)];
if fitEdge
	A3 = photoElectric(E)./photoElectric(E0) .* HKappa(E, EKappa, epsilonKappa);
	A = [A, A3];
end
if modelPairProduction
	A = [A, pairProduction(E)];
end

% weight errors with inverse density of given data points
deltaE = diff(E);
weightsDensity = zeros(size(E));
weightsDensity(1:end-1) = deltaE/2;
weightsDensity(2:end) = weightsDensity(2:end) + deltaE/2;
weightsDensity(1) = 2*weightsDensity(1);
weightsDensity(end) = 2*weightsDensity(end);

% weight errors with energy (since detector intensity is mostly proportional to photon energy)
weightsEnergy = E;

% weight errors with spectrum
% weightsSpectrum = ;

% extra weight for specific energy levels
weightsSpecificE = ones(size(E));
weightsSpecificE(E == E0) = 100 * weightsSpecificE(E == E0);
weightsSpecificE(E == EMeV) = 10 * weightsSpecificE(E == EMeV);

% accumulate weights and create scaling matrix
weights = weightsDensity .* weightsSpecificE;
% weights = weightsDensity .* weightsEnergy .* weightsSpectrum .* weightsSpecificE;

% apply weights to problem
b = diag(weights) * b;
A = diag(weights) * A;

% optimize and store material parameters
x = A \ b;
Phi = x(1);
Theta = x(2);
if fitEdge
	Kappa = x(3);
else
	Kappa = 0;
end
if modelPairProduction
	if fitEdge
		Pi = x(4);
	else
		Pi = x(3);
	end
else
	Pi = 0;
end


%% Plot Comparison between Actual and Approximated Attenuation Profile

if interactive
	% compute actual and approximated attenuation at E0
	muApprox = Phi * photoElectric(E)./photoElectric(E0) + Theta * kleinNishina(E)./kleinNishina(E0) + Kappa * photoElectric(E)./photoElectric(E0) .* HKappa(E, EKappa, epsilonKappa) + Pi * pairProduction(E);

	% plot
	figure('Name', ['Attenuation of ', materialName], 'NumberTitle', 'off');
	zoom on;
	loglog(E, mu, 'b-');
	hold on;
	loglog(E, muApprox, 'r-');
	loglog(E0, muE0, 'bx', 'MarkerSize', 10);
	hold off;
	xlabel('Photon energy (keV)');
	ylabel('Attenuation (1/cm)');
	legend('mu', 'mu_{approx}');
	grid on;
end
