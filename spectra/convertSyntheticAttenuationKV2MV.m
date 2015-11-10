function imgAttMeV = convertSyntheticAttenuationKV2MV( imgAttKeV, spectrumKeV, spectrumMeV, useMonochromatic )
% convert attenuation in KeV image to the corresponding values in MeV image
% input:
%       imgAttKeV   - attenuation in KeV image
%       spectrumKeV - X-ray spectrum infomation of KeV
%       spectrumMeV - X-ray spectrum infomation of MeV
% output:
%       imgAttMeV   - attenuation in MeV image
%
% Meng Wu at Stanford
% 2012.10

if nargin < 4
    useMonochromatic = false;
end

% Compute attenuations lookup table for two sprectrums
materialNames = {'Air_Dry_near_sea_level', ...
    'Adipose_Tissue_ICRU-44', ...
    'Water_liquid', ...
    'Muscle_Skeletal_ICRU-44', ...
    'Tissue_Soft_ICRU-44', ...
    'B-100_Bone-Equivalent_Plastic', ...
    'Bone_Cortical_ICRU-44',...
    'Aluminum',...
    'Titanium',...
    'Iron', ...
    'Copper'};

HighestAttMaterial = 5;

muKeV = zeros(HighestAttMaterial, 1);
muMeV = zeros(HighestAttMaterial, 1);

for m = 1 : HighestAttMaterial

    
    if useMonochromatic %monochromatic
        muKeV(m) = materialAttenuation( spectrumKeV.energyAverage, materialNames{m});
        muMeV(m) = materialAttenuation( spectrumMeV.energyAverage, materialNames{m});
    else
        [~, muKeV(m)] = materialAttenuation( spectrumKeV.energyBinLabels, ...
            materialNames{m}, spectrumKeV.photonsPerEnergyBin);
        [~, muMeV(m)] = materialAttenuation( spectrumMeV.energyBinLabels, ...
            materialNames{m}, spectrumMeV.photonsPerEnergyBin);
    end
    
end

% make sure no pixel is outside the lookup table
imgAttKeV( imgAttKeV < muKeV(1) ) = muKeV(1);
%imgAttKeV( imgAttKeV > muKeV(end)) =  muKeV(end);

% attenuation conversion using lookup table
imgAttMeV = imgAttKeV;
imgAttMeV( imgAttKeV < muKeV(end) ) = interp1( muKeV, muMeV, imgAttKeV( imgAttKeV < muKeV(end)));

