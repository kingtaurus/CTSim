function [...
	imgEYPhi, imgEYTheta, imgEYKappa, imgEYPi, ...
	imgM, imgN, imgP, imgQ] ...
	= reconImpactIntermediates(...
	SAD, ADD, sizeDet, spacingDet, oversamplingDetector, overSamplingDepth, betas, ...
	energyBinLabels, photonsPerEnergyBin, sinoPhotonsNoisy, sinoMapWork, ...
	E0, imgPhi, imgTheta, imgKappa, imgPi, imgPhiPrime, imgThetaPrime, imgKappaPrime, imgPiPrime, ...
	modelEdge, EKappa, epsilonKappa, modelPairProduction, ...
	spacingRecon, ...
	outputDir)
% computes all the intermediate variables used in an iteration of
% Bruno De Man's IMPACT algorithm (Trans. Med. Imag. 20(10), pp. 999-1008, 2001)
%
% Copyright (c) 2012 by Andreas Keil, Stanford University.



%% Validate/Parse Arguments

validateattributes(SAD, {'numeric'}, {'scalar', 'real', 'positive'});
validateattributes(ADD, {'numeric'}, {'scalar', 'real', 'nonnegative'});
validateattributes(sizeDet, {'numeric'}, {'scalar', 'integer', 'positive'});
validateattributes(spacingDet, {'numeric'}, {'scalar', 'real', 'positive'});
validateattributes(oversamplingDetector, {'numeric'}, {'scalar', 'integer', 'positive'});
validateattributes(overSamplingDepth, {'numeric'}, {'scalar', 'real', 'positive'});
validateattributes(betas, {'numeric'}, {'vector', 'real'});
sizeSino = [sizeDet length(betas)];
validateattributes(energyBinLabels, {'numeric'}, {'vector', 'real', 'nonnegative'});
validateattributes(photonsPerEnergyBin, {'numeric'}, {'size', [length(energyBinLabels) 1], 'real', 'nonnegative'});
validateattributes(sinoPhotonsNoisy, {'numeric'}, {'size', sizeSino, 'real', 'nonnegative'});
if isempty(sinoMapWork)
	sinoMapWork = true(size(sinoPhotonsNoisy));
else
	validateattributes(sinoMapWork, {'logical'}, {'size', sizeSino});
end
validateattributes(E0, {'numeric'}, {'scalar', 'real', 'positive'});
validateattributes(imgPhi, {'numeric'}, {'2d', 'real'});
sizeRecon = size(imgPhi);
validateattributes(imgTheta, {'numeric'}, {'size', sizeRecon, 'real'});
validateattributes(imgPhiPrime, {'numeric'}, {'size', sizeRecon, 'real'});
validateattributes(imgThetaPrime, {'numeric'}, {'size', sizeRecon, 'real'});
validateattributes(modelEdge, {'logical'}, {'scalar'});
if modelEdge
	validateattributes(EKappa, {'numeric'}, {'scalar', 'real', 'positive'});
	validateattributes(epsilonKappa, {'numeric'}, {'scalar', 'real', 'positive'});
	validateattributes(imgKappa, {'numeric'}, {'size', sizeRecon, 'real'});
	validateattributes(imgKappaPrime, {'numeric'}, {'size', sizeRecon, 'real'});
end
if modelPairProduction
	validateattributes(imgPi, {'numeric'}, {'size', sizeRecon, 'real'});
	validateattributes(imgPiPrime, {'numeric'}, {'size', sizeRecon, 'real'});
end
validateattributes(spacingRecon, {'numeric'}, {'size', [1 2], 'real', 'positive'});
validateattributes(outputDir, {'char'}, {});


%% Forward Project Coefficient Volumes

sinoPhi = forwardProjectSetOfViews(imgPhi, spacingRecon, SAD, ADD, betas, sizeDet, spacingDet, oversamplingDetector, overSamplingDepth, sinoMapWork);
sinoTheta = forwardProjectSetOfViews(imgTheta, spacingRecon, SAD, ADD, betas, sizeDet, spacingDet, oversamplingDetector, overSamplingDepth, sinoMapWork);
if modelEdge, sinoKappa = forwardProjectSetOfViews(imgKappa, spacingRecon, SAD, ADD, betas, sizeDet, spacingDet, oversamplingDetector, overSamplingDepth, sinoMapWork); else sinoKappa = 0; end
if modelPairProduction, sinoPi = forwardProjectSetOfViews(imgPi, spacingRecon, SAD, ADD, betas, sizeDet, spacingDet, oversamplingDetector, overSamplingDepth, sinoMapWork); else sinoPi = 0; end


%% Accumulate Intermediate Variables over Energy Bins

yhat = zeros(size(sinoPhi));
YPhi = zeros(size(sinoPhi));
YTheta = zeros(size(sinoPhi));
YPhiPhi = zeros(size(sinoPhi));
YPhiTheta = zeros(size(sinoPhi));
YThetaTheta = zeros(size(sinoPhi));
if modelEdge
	YKappa = zeros(size(sinoPhi));
	YPhiKappa = zeros(size(sinoPhi));
	YThetaKappa = zeros(size(sinoPhi));
	YKappaKappa = zeros(size(sinoPhi));
end
if modelPairProduction
	YPi = zeros(size(sinoPhi));
	YPhiPi = zeros(size(sinoPhi));
	YThetaPi = zeros(size(sinoPhi));
	YPiPi = zeros(size(sinoPhi));
end
if modelEdge && modelPairProduction
	YKappaPi = zeros(size(sinoPhi));
end
for k = 1:length(energyBinLabels)
	Phi_k = photoElectric(energyBinLabels(k))./photoElectric(E0);
	Theta_k = kleinNishina(energyBinLabels(k))./kleinNishina(E0);
	if modelEdge, Kappa_k = photoElectric(energyBinLabels(k))./photoElectric(E0) * HKappa(energyBinLabels(k), EKappa, epsilonKappa); else Kappa_k = 0; end
	if modelPairProduction, Pi_k = pairProduction(energyBinLabels(k)); else Pi_k = 0; end
	yhat_k = photonsPerEnergyBin(k)*exp(-sinoPhi*Phi_k-sinoTheta*Theta_k-sinoKappa*Kappa_k-sinoPi*Pi_k);
	yhat = yhat + yhat_k;
	YPhi = YPhi + Phi_k*yhat_k;
	YTheta = YTheta + Theta_k*yhat_k;
	YPhiPhi = YPhiPhi + Phi_k^2*yhat_k;
	YPhiTheta = YPhiTheta + Phi_k*Theta_k*yhat_k;
	YThetaTheta = YThetaTheta + Theta_k^2*yhat_k;
	if modelEdge
		YKappa = YKappa + Kappa_k*yhat_k;
		YPhiKappa = YPhiKappa + Phi_k*Kappa_k*yhat_k;
		YThetaKappa = YThetaKappa + Theta_k*Kappa_k*yhat_k;
		YKappaKappa = YKappaKappa + Kappa_k^2*yhat_k;
	end
	if modelPairProduction
		YPi = YPi + Pi_k*yhat_k;
		YPhiPi = YPhiPi + Phi_k*Pi_k*yhat_k;
		YThetaPi = YThetaPi + Theta_k*Pi_k*yhat_k;
		YPiPi = YPiPi + Pi_k^2*yhat_k;
	end
	if modelEdge && modelPairProduction
		YKappaPi = YKappaPi + Kappa_l*Pi_k*yhat_k;
	end
end
clear sinoPhi sinoTheta sinoKappa sinoPi yhat_k Phi_k Theta_k Kappa_k Pi_k;


%% Compute Relative Pixel Errors

sinoE = 1 - sinoPhotonsNoisy./yhat;
sinoE(yhat == 0) = 0;


%% Forward Project Coefficient Derivates

sinoPhiPrime = forwardProjectSetOfViews(imgPhiPrime, spacingRecon, SAD, ADD, betas, sizeDet, spacingDet, oversamplingDetector, 2, sinoMapWork);
sinoThetaPrime = forwardProjectSetOfViews(imgThetaPrime, spacingRecon, SAD, ADD, betas, sizeDet, spacingDet, oversamplingDetector, 2, sinoMapWork);
if modelEdge, sinoKappaPrime = forwardProjectSetOfViews(imgKappaPrime, spacingRecon, SAD, ADD, betas, sizeDet, spacingDet, oversamplingDetector, 2, sinoMapWork); end
if modelPairProduction, sinoPiPrime = forwardProjectSetOfViews(imgPiPrime, spacingRecon, SAD, ADD, betas, sizeDet, spacingDet, oversamplingDetector, 2, sinoMapWork); end


%% Compute Auxiliary Variables for the Denominator's Intermediate Variables

%TODO: maybe save memory by discarding yi, yi^hat, Y^PhiPhi, Y^ThetaTheta asap
M = sinoPhiPrime .* (YPhiPhi.*sinoE + sinoPhotonsNoisy.*YPhi.^2./yhat.^2) ...
	+ sinoThetaPrime .* (YPhiTheta.*sinoE + sinoPhotonsNoisy.*YPhi.*YTheta./yhat.^2);
N = sinoPhiPrime .* (YPhiTheta.*sinoE + sinoPhotonsNoisy.*YPhi.*YTheta./yhat.^2) ...
	+ sinoThetaPrime .* (YThetaTheta.*sinoE + sinoPhotonsNoisy.*YTheta.^2./yhat.^2);
if modelEdge
	M = M + sinoKappaPrime .* (YPhiKappa.*sinoE + sinoPhotonsNoisy.*YPhi.*YKappa./yhat.^2);
	N = N + sinoKappaPrime .* (YThetaKappa.*sinoE + sinoPhotonsNoisy.*YTheta.*YKappa./yhat.^2);
	P = sinoPhiPrime .* (YPhiKappa.*sinoE + sinoPhotonsNoisy.*YPhi.*YKappa./yhat.^2) ...
		+ sinoThetaPrime .* (YThetaKappa.*sinoE + sinoPhotonsNoisy.*YTheta.*YKappa./yhat.^2) ...
		+ sinoKappaPrime .* (YKappaKappa.*sinoE + sinoPhotonsNoisy.*YKappa.^2./yhat.^2);
end
if modelPairProduction
	M = M + sinoPiPrime .* (YPhiPi.*sinoE + sinoPhotonsNoisy.*YPhi.*YPi./yhat.^2);
	N = N + sinoPiPrime .* (YThetaPi.*sinoE + sinoPhotonsNoisy.*YTheta.*YPi./yhat.^2);
	Q = sinoPhiPrime .* (YPhiPi.*sinoE + sinoPhotonsNoisy.*YPhi.*YPi./yhat.^2) ...
		+ sinoThetaPrime .* (YThetaPi.*sinoE + sinoPhotonsNoisy.*YTheta.*YPi./yhat.^2) ...
		+ sinoPiPrime .* (YPiPi.*sinoE + sinoPhotonsNoisy.*YPi.^2./yhat.^2);
end
if modelEdge && modelPairProduction
	P = P + sinoPiPrime .* (YKappaPi.*sinoE + sinoPhotonsNoisy.*YKappa.*YPi./yhat.^2);
	Q = Q + sinoKappaPrime .* (YKappaPi.*sinoE + sinoPhotonsNoisy.*YKappa.*YPi./yhat.^2);
end
clear sinoPhiPrime sinoThetaPrime sinoKappaPrime sinoPiPrime YPhiPhi YPhiTheta YPhiKappa YPhiPi YThetaTheta YThetaKappa YThetaPi YPiPi yhat;


%% Copmute Back Projections for the Intermediate Variables that Are to Be Returned

imgEYPhi = backProjectPlainSetOfViews(sinoE.*YPhi, betas, spacingDet, SAD, ADD, sizeRecon, spacingRecon, false, sinoMapWork);
imgEYTheta = backProjectPlainSetOfViews(sinoE.*YTheta, betas, spacingDet, SAD, ADD, sizeRecon, spacingRecon, false, sinoMapWork);
if modelEdge, imgEYKappa = backProjectPlainSetOfViews(sinoE.*YKappa, betas, spacingDet, SAD, ADD, sizeRecon, spacingRecon, false, sinoMapWork); else imgEYKappa = []; end
if modelPairProduction, imgEYPi = backProjectPlainSetOfViews(sinoE.*YPi, betas, spacingDet, SAD, ADD, sizeRecon, spacingRecon, false, sinoMapWork); else imgEYPi = []; end
imgM = backProjectPlainSetOfViews(M, betas, spacingDet, SAD, ADD, sizeRecon, spacingRecon, false, sinoMapWork);
imgN = backProjectPlainSetOfViews(N, betas, spacingDet, SAD, ADD, sizeRecon, spacingRecon, false, sinoMapWork);
if modelEdge, imgP = backProjectPlainSetOfViews(P, betas, spacingDet, SAD, ADD, sizeRecon, spacingRecon, false, sinoMapWork); else imgP = []; end
if modelPairProduction, imgQ = backProjectPlainSetOfViews(Q, betas, spacingDet, SAD, ADD, sizeRecon, spacingRecon, false, sinoMapWork); else imgQ = []; end
