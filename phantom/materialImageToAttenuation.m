function imgAtt = materialImageToAttenuation(phan, spectrum)
% Generate Attenuation Volume


imgMat = phan.phantomMaterials;
spectrumEnergies = spectrum.energyAverage;
photonsPerEnergyBin = spectrum.photonsTotal;
materialIndexes = phan.materialIndexes;
materials = phan.materialNames;

imgAtt = zeros(size(imgMat));

for m = 1:length(materials)
    
    [mu, muEff] =  materialAttenuation( spectrumEnergies, ...
        materials{m}, photonsPerEnergyBin);
    
    imgAtt(imgMat == materialIndexes{m}) = muEff;
    
end
if phan.useDensity
    imgAtt = imgAtt .* phan.phantomDensities;
end
