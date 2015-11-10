function polyCoeffs = beamHardeningCorrectionPolynomialCoefficients( material, ...
    energyBinLabels, photonsPerEnergyBin, energyEff, order,  maximumThickness )

if nargin < 6
    maximumThickness = 50;
end

photonsTotal = sum( photonsPerEnergyBin );

[attenuationPoly, attenuationMono] = materialAttenuation( energyBinLabels, material, photonsPerEnergyBin, energyEff);

% initialize vectors
intensityPolyRef = zeros(256, 1);
intensityMonoRef = zeros(256, 1);


intersectionLengths = linspace(0,maximumThickness,256);

for i = 1:256
    l = intersectionLengths(i);
    intensityPolyRef(i) = sum(photonsPerEnergyBin.*exp(-l*attenuationPoly));
    intensityMonoRef(i) = photonsTotal*exp(-l*attenuationMono);
end

lineIntegralPoly = -log(intensityPolyRef / photonsTotal);
lineIntegralMono = -log(intensityMonoRef / photonsTotal);

% least squares computing polynomial coefficients
A = ones( 256, order  );
for i = 1 : order
    A(:,i) = lineIntegralPoly.^i;
end

polyCoeffs = ( A' * A ) \ ( A' * lineIntegralMono );


end