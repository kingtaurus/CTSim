function sinoOut = beamHardeningWarterCorrection(sinoIn, spectrum, multiChannels, polyOrders, maximumThickness )
% function sinoOut = beamHardeningWarterCorrection(sinoIn, spectrum, multiChannels, polyOrders )
%
% multi-channal polynomial water based beam hardening correction
%
% Meng Wu at Stanford University
% 2013


if nargin < 3
    multiChannels = false;
end

if nargin < 4
    polyOrders = 3;
end

if nargin < 5
    maximumThickness = 50;
end

% scale the spectrum if energy integrating detector
if spectrum.energyIntegrating
    photonsPerEnergyBin = spectrum.DQE * spectrum.photonsPerEnergyBinOriginal .* spectrum.energyBinLabels * spectrum.detectorGain ;
else
    photonsPerEnergyBin = spectrum.DQE * spectrum.photonsPerEnergyBinOriginal * spectrum.detectorGain  ;
end

if ~ multiChannels

    if spectrum.useBowtie
        photonsPerEnergyBinChannel = computeResultingPhotons( photonsPerEnergyBin, ...
            spectrum.energyBinLabels, spectrum.bowtieMaterial, mean( spectrum.bowtieThickness(:)) );
    end
    
    polyCoeffs = beamHardeningCorrectionPolynomialCoefficients(  'Water_Liquid', ...
        spectrum.energyBinLabels, photonsPerEnergyBinChannel, spectrum.energyAverage, polyOrders, maximumThickness );
    
    sinoOut = sinoIn * polyCoeffs(1);
    for i = 2 : polyOrders
        sinoOut = sinoOut + sinoIn.^i * polyCoeffs(i);
    end
    
else
    
    numChannels = size( spectrum.flatFieldRatio, 2);
    %polyCoeffsTable = zeros( polyOrders, numChannels );
    sinoOut = sinoIn;
    
    for ic = 1 : numChannels
        photonsPerEnergyBinChannel = computeResultingPhotons( photonsPerEnergyBin, ...
            spectrum.energyBinLabels, spectrum.bowtieMaterial, spectrum.bowtieThickness(1,ic) );
        
        polyCoeffs = beamHardeningCorrectionPolynomialCoefficients(  'Water_Liquid', ...
            spectrum.energyBinLabels, photonsPerEnergyBinChannel, spectrum.energyAverage, polyOrders );
        
        if length( size(sinoIn) ) == 3
            sinoOut(:,ic,:) = sinoIn(:,ic,:)  * polyCoeffs(1);
            for i = 2 : polyOrders
                sinoOut(:,ic,:) = sinoOut(:,ic,:)  + sinoIn(:,ic,:) .^i * polyCoeffs(i);
            end
        else
            sinoOut(ic,:) = sinoIn(ic,:)  * polyCoeffs(1);
            for i = 2 : polyOrders
                sinoOut(ic,:)  = sinoOut(ic,:)  + sinoIn(ic,:) .^i * polyCoeffs(i);
            end
        end
        
    end
    
end


end