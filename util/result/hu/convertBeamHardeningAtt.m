function out = convertBeamHardeningAtt( in,spectrum, mode )
% Compute Beam Hardening Correction in sinogram space using Joseph's method.
%   input:
%       in       - sinogram or image in attenuation
%       spectrum - spectrum infomation
%       mode     - 0. assuming ploychromatic image (default)
%                  1. assuming monochromatic image
%   output:
%       sinoOut  - correction to add to original sinogram

if nargin < 3
    mode = 0;
end

% Compute attenuations lookup table for two sprectrums
materialsDir = 'physicsdata/materials/';

materials = {'Air_Dry_near_sea_level', ...
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

[E, mu] = readMaterialAttenuations(materials, materialsDir);

noMaterials = length(materials);

muPoly = zeros(noMaterials, 1);
muMono = zeros(noMaterials, 1);

for m = 1 : noMaterials
    
    muMono(m) = spectralAttenuation(spectrum.energyAverage, 1, E{m}, mu{m});
    
    muPoly(m) = spectralAttenuation(spectrum.energyBinLabels, ...
        spectrum.photonsPerEnergyBin, E{m}, mu{m});
    
end


% make sure no pixel is outside the lookup table
if mode == 0
    
    in( in < muPoly(2) * 0.3 ) = muPoly(1);
    in( in > muPoly(noMaterials)) =  muPoly(noMaterials);
    out = interp1( muPoly, muMono, in);
    
else
    
    in( in < muMono(2) * 0.3 ) = muMono(1);
    in( in > muMono(noMaterials)) =  muMono(noMaterials);
    out = interp1( muMono, muPoly, in);
    
end

end