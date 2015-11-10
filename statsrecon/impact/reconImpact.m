function imgReconAtt = reconImpact(...
	SAD, ADD, sizeDet, spacingDet, oversamplingDetector, overSamplingDepth, betas, ...
	energyBinLabels, photonsPerEnergyBin, sinoPhotonsNoisy, sinoMapWork, ...
	E0, deManMu, deManPhi, deManTheta, deManKappa, deManPi, ...
	modelEdge, EKappa, epsilonKappa, modelPairProduction, ...
	imgReconAtt, spacingRecon, regularizationWeight, regularizationParam, weightUpdate, ...
	saveIntermediateSteps, outputDir)
% performs an iteration of Bruno De Man's IMPACT algorithm
% (Trans. Med. Imag. 20(10), pp. 999-1008, 2001)
%
% Copyright (c) 2012 by Andreas Keil, Stanford University.



%% Validate/Parse Arguments

validateattributes(SAD, {'cell'}, {'vector'});
noGeometries = length(SAD);
sizeSino = cell(noGeometries, 1);
for g = 1:noGeometries
	validateattributes(SAD{g}, {'numeric'}, {'scalar', 'real', 'positive'});
	validateattributes(ADD{g}, {'numeric'}, {'scalar', 'real', 'nonnegative'});
	validateattributes(sizeDet{g}, {'numeric'}, {'scalar', 'integer', 'positive'});
	validateattributes(spacingDet{g}, {'numeric'}, {'scalar', 'real', 'positive'});
	validateattributes(betas{g}, {'numeric'}, {'vector', 'real'});
	sizeSino{g} = [sizeDet{g} length(betas{g})];
	validateattributes(energyBinLabels{g}, {'numeric'}, {'vector', 'real', 'nonnegative'});
	validateattributes(photonsPerEnergyBin{g}, {'numeric'}, {'size', [length(energyBinLabels{g}) 1], 'real', 'nonnegative'});
	validateattributes(sinoPhotonsNoisy{g}, {'numeric'}, {'size', sizeSino{g}, 'real', 'nonnegative'});
	if isempty(sinoMapWork{g})
		sinoMapWork{g} = true(size(sinoPhotonsNoisy{g}));
	else
		validateattributes(sinoMapWork{g}, {'logical'}, {'size', sizeSino{g}});
	end
end
validateattributes(oversamplingDetector, {'numeric'}, {'scalar', 'integer', 'positive'});
validateattributes(overSamplingDepth, {'numeric'}, {'scalar', 'real', 'positive'});
validateattributes(E0, {'numeric'}, {'scalar', 'real', 'positive'});
validateattributes(deManMu, {'numeric'}, {'vector', 'real', 'nonnegative'});
noDeManSamplePoints = length(deManMu);
validateattributes(deManPhi, {'numeric'}, {'size', [noDeManSamplePoints 1], 'real', 'nonnegative'});
validateattributes(deManTheta, {'numeric'}, {'size', [noDeManSamplePoints 1], 'real', 'nonnegative'});
validateattributes(modelEdge, {'logical'}, {'scalar'});
if modelEdge
	validateattributes(deManKappa, {'numeric'}, {'size', [noDeManSamplePoints 1], 'real', 'nonnegative'});
	validateattributes(EKappa, {'numeric'}, {'scalar', 'real', 'positive'});
	validateattributes(epsilonKappa, {'numeric'}, {'scalar', 'real', 'positive'});
end
if modelPairProduction
	validateattributes(deManPi, {'numeric'}, {'size', [noDeManSamplePoints 1], 'real', 'nonnegative'});
end
validateattributes(imgReconAtt, {'numeric'}, {'2d', 'real'});
validateattributes(spacingRecon, {'numeric'}, {'size', [1 2], 'real', 'positive'});
validateattributes(saveIntermediateSteps, {'logical'}, {'scalar'});
validateattributes(outputDir, {'char'}, {});


%% Compute Coefficient Volumes

imgPhi = interp1(deManMu, deManPhi, imgReconAtt, 'linear', 'extrap');
if saveIntermediateSteps, writeImageFormats(imgPhi, spacingRecon, [min(imgPhi(:)), max(imgPhi(:))], [outputDir 'ReconImpact_imgPhi']); end
imgTheta = interp1(deManMu, deManTheta, imgReconAtt, 'linear', 'extrap');
if saveIntermediateSteps, writeImageFormats(imgTheta, spacingRecon, [min(imgTheta(:)), max(imgTheta(:))], [outputDir 'ReconImpact_imgTheta']); end
if modelEdge
	imgKappa = interp1(deManMu, deManKappa, imgReconAtt, 'linear', 'extrap');
	if saveIntermediateSteps, writeImageFormats(imgKappa, spacingRecon, [min(imgKappa(:)), max(imgKappa(:))], [outputDir 'ReconImpact_imgKappa']); end
else
	imgKappa = [];
end
if modelPairProduction
	imgPi = interp1(deManMu, deManPi, imgReconAtt, 'linear', 'extrap');
	if saveIntermediateSteps, writeImageFormats(imgPi, spacingRecon, [min(imgPi(:)), max(imgPi(:))], [outputDir 'ReconImpact_imgPi']); end
else
	imgPi = [];
end


%% Compute Coefficient Derivative Volumes

% compute derivatives for all sections of the piecewise linear coefficient
% functions and make sure that there is a somewhat smooth transition
% between sections by defining two sampling points for each given mu value
deManMuForPrime = [deManMu(:)'/(1+10*eps); deManMu(:)'*(1+10*eps)]; deManMuForPrime = deManMuForPrime(2:end-1)';
deManPhiPrime = diff(deManPhi)./diff(deManMu); deManPhiPrime = [deManPhiPrime(:)'; deManPhiPrime(:)']; deManPhiPrime = deManPhiPrime(:);
deManThetaPrime = diff(deManTheta)./diff(deManMu); deManThetaPrime = [deManThetaPrime(:)'; deManThetaPrime(:)']; deManThetaPrime = deManThetaPrime(:);
if modelEdge, deManKappaPrime = diff(deManKappa)./diff(deManMu); deManKappaPrime = [deManKappaPrime(:)'; deManKappaPrime(:)']; deManKappaPrime = deManKappaPrime(:); end
if modelPairProduction, deManPiPrime = diff(deManPi)./diff(deManMu); deManPiPrime = [deManPiPrime(:)'; deManPiPrime(:)']; deManPiPrime = deManPiPrime(:); end

% compute coefficient derivates at the mu values of the current image
imgPhiPrime = interp1(deManMuForPrime, deManPhiPrime, imgReconAtt, 'linear', 'extrap');
imgThetaPrime = interp1(deManMuForPrime, deManThetaPrime, imgReconAtt, 'linear', 'extrap');
if modelEdge, imgKappaPrime = interp1(deManMuForPrime, deManKappaPrime, imgReconAtt, 'linear', 'extrap'); else imgKappaPrime = []; end
if modelPairProduction, imgPiPrime = interp1(deManMuForPrime, deManPiPrime, imgReconAtt, 'linear', 'extrap'); else imgPiPrime = []; end

clear deManMuForPrime deManPhiPrime deManThetaPrime deManKappaPrime deManPiPrime;


%% Compute all Intermediate Variables of Numerator and Denominator

% compute variables per geometry
imgEYPhi{g} = cell(noGeometries, 1);
imgEYTheta{g} = cell(noGeometries, 1);
if modelEdge, imgEYKappa{g} = cell(noGeometries, 1); else imgEYKappa{g} = []; end
if modelPairProduction, imgEYPi{g} = cell(noGeometries, 1); else imgEYPi{g} = []; end
imgM{g} = cell(noGeometries, 1);
imgN{g} = cell(noGeometries, 1);
if modelEdge, imgP{g} = cell(noGeometries, 1); else imgP{g} = []; end
if modelPairProduction, imgQ{g} = cell(noGeometries, 1); else imgQ{g} = []; end
for g = 1:noGeometries
	[...
		imgEYPhi{g}, imgEYTheta{g}, imgEYKappa{g}, imgEYPi{g}, ...
		imgM{g}, imgN{g}, imgP{g}, imgQ{g}] ...
		= reconImpactIntermediates(...
			SAD{g}, ADD{g}, sizeDet{g}, spacingDet{g}, oversamplingDetector, overSamplingDepth, betas{g}, ...
			energyBinLabels{g}, photonsPerEnergyBin{g}, sinoPhotonsNoisy{g}, sinoMapWork{g}, ...
			E0, imgPhi, imgTheta, imgKappa, imgPi, imgPhiPrime, imgThetaPrime, imgKappaPrime, imgPiPrime, ...
			modelEdge, EKappa, epsilonKappa, ...
			modelPairProduction, ...
			spacingRecon, ...
			outputDir...
		);
end

% accumulate variables
sizeRecon = size(imgReconAtt);
imgEYPhiSum = zeros(sizeRecon);
imgEYThetaSum = zeros(sizeRecon);
if modelEdge, imgEYKappaSum = zeros(sizeRecon); else imgEYKappaSum = []; end
if modelPairProduction, imgEYPiSum = zeros(sizeRecon); else imgEYPiSum = []; end
imgMSum = zeros(sizeRecon);
imgNSum = zeros(sizeRecon);
if modelEdge, imgPSum = zeros(sizeRecon); else imgPSum = []; end
if modelPairProduction, imgQSum = zeros(sizeRecon); else imgQSum = []; end
for g = 1:noGeometries
	imgEYPhiSum = imgEYPhiSum + imgEYPhi{g};
	imgEYThetaSum = imgEYThetaSum + imgEYTheta{g};
	if modelEdge, imgEYKappaSum = imgEYKappaSum + imgEYKappa{g}; end
	if modelPairProduction, imgEYPiSum = imgEYPiSum + imgEYPi{g}; end
	imgMSum = imgMSum + imgM{g};
	imgNSum = imgNSum + imgN{g};
	if modelEdge, imgPSum = imgPSum + imgP{g}; end
	if modelPairProduction, imgQSum = imgQSum + imgQ{g}; end
end


%% Assemble Update

deltaAttNumerator = imgPhiPrime.*imgEYPhiSum + imgThetaPrime.*imgEYThetaSum;
deltaAttDenominator = imgPhiPrime.*imgMSum + imgThetaPrime.*imgNSum;
if modelEdge
	deltaAttNumerator = deltaAttNumerator + imgKappaPrime.*imgEYKappaSum;
	deltaAttDenominator = deltaAttDenominator + imgKappaPrime.*imgPSum;
end
if modelPairProduction
	deltaAttNumerator = deltaAttNumerator + imgPiPrime.*imgEYPiSum;
	deltaAttDenominator = deltaAttDenominator + imgPiPrime.*imgQSum;
end
if regularizationWeight > 0
	deltaMu = imgReconAtt(2:end, :) - imgReconAtt(1:end-1, :);
	deltaAttNumerator(1:end-1, :) = deltaAttNumerator(1:end-1, :) + ...
		regularizationWeight * ...
		dHuber(deltaMu, regularizationParam);
	deltaMu = imgReconAtt(1:end-1, :) - imgReconAtt(2:end, :);
	deltaAttNumerator(2:end, :) = deltaAttNumerator(2:end, :) + ...
		regularizationWeight * ...
		dHuber(deltaMu, regularizationParam);
	deltaMu = imgReconAtt(:, 2:end) - imgReconAtt(:, 1:end-1);
	deltaAttNumerator(:, 1:end-1) = deltaAttNumerator(:, 1:end-1) + ...
		regularizationWeight * ...
		dHuber(deltaMu, regularizationParam);
	deltaMu = imgReconAtt(:, 1:end-1) - imgReconAtt(:, 2:end);
	deltaAttNumerator(:, 2:end) = deltaAttNumerator(:, 2:end) + ...
		regularizationWeight * ...
		dHuber(deltaMu, regularizationParam);
end
deltaAtt = deltaAttNumerator./deltaAttDenominator;


%% Prevent NaNs for Image Pixles that Did not Get Any Contribution from Any Projection (Which Can Happen for Subsets)

deltaAtt(deltaAttDenominator == 0) = 0;


%% Prevent Huge Updates

signDelta = sign(deltaAtt);
absDelta = abs(deltaAtt);
cap1 = 200 * median(absDelta(:));
cap2 = 0.1 * abs(imgReconAtt);
mapCap1 = absDelta > cap1;
mapCap2 = absDelta > cap2;
mapCapBoth = mapCap1 & mapCap2;
deltaAtt(mapCapBoth) = (cap1 + cap2(mapCapBoth))/2 .* signDelta(mapCapBoth);

% signDelta = sign(deltaAtt);
% absDeltaRel = abs(deltaAtt ./ (imgReconAtt + 1e-6));
% cap = 10 * median(absDeltaRel(:));
% mapCap = absDeltaRel > cap;
% deltaAtt(mapCap) = cap * signDelta(mapCap) .* imgReconAtt(mapCap);


%% Final Check that Nothing Went Wrong

bigUpdates = logical(abs(deltaAtt) > 0.25);
if any(bigUpdates(:))
    warning('\nUpdate is really big.\n'); %#ok<WNTAG>
    figure; spy(bigUpdates);
end
if any(isnan(deltaAtt(:)))
	warning('\nISNAN!\n'); %#ok<WNTAG>
end


%% Apply the Update

imgReconAtt = imgReconAtt + weightUpdate*deltaAtt;
