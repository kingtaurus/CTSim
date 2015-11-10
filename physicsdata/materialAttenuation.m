function [mu, muEff] = materialAttenuation( energy, material, photonsPerEnergy, energyEff)
% Return linear attenuation coefficeints of materails at given x-ray energy
% Inputs: 
%       energy - X-ray photon energies
%       material - material name, see below for defined material names
%       photonsPerEnergy    - use for computing effective energy
%       energyEff           - user defined effective energy
% Outputs:
%       mu  - attenuation coeffient at each photon energy (cm^-1)
%       muEff  - effective attenuation coeffient (cm^-1)
%
% NIST materials:
%       {'Air_Dry_near_sea_level', 'Adipose_Tissue_ICRU-44', 'Water_Liquid', ...
%       'Tissue_Soft_ICRU-44', 'Muscle_Skeletal_ICRU-44', 'B-100_Bone-Equivalent_Plastic', ...
%       'Bone_Cortical_ICRU-44', 'Aluminum', 'Titanium', 'Iron', 'Copper', 'AK-Amalgam', 'Gold'}
% Self defined materials:
%       {'adipose', 'blood', 'bone_compact', 'bone_cortical', 'brain', 'lung', ...
%        'muscle_skeletal', 'muscle_striated', 'skin', 'soft_tissue', 'water', 'air' ...
%        'CWO', 'Acrylic', 'PMP', 'Delrin', 'Teflon', 'Polystryrene', 'Bone20', 'Bone50' ...,
%        'LDPE' };
% Elements with alphabetical symbles
%
% Meng Wu at Stanford University
% 2013 - 2014

if findNISTMaterial( material )
    [E, mu_E] = readMaterialAttenuation(material, 'physicsdata/materials/');
    mu = spectralAttenuations(energy, E, mu_E);
else
    [ztable, density ] = findCompositionMaterial( material );
    mu = density * XrayMu(ztable, energy) ;
end

if length( energy ) == 1
    muEff = mu;

elseif nargin == 2
    muEff = median( mu );
elseif nargin == 3
    energyEff = sum( energy(:) .* photonsPerEnergy(:) )/ sum( photonsPerEnergy );
    muEff = interp1geom(energy, mu, energyEff);
else
    muEff = interp1geom(energy, mu, energyEff);
end


end

