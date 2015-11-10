function [massAtt, massAttEff] = materialMassAttenuation( energy, material, photonsPerEnergy, energyEff)
% Return mass attenuation coefficeints of materails at given x-ray energy 
% ( cm^2 / g )
% Inputs: 
%       energy - X-ray photon energies
%       material - material name, see below for defined material names
%       photonsPerEnergy    - use for computing effective energy
%       energyEff           - user defined effective energy
% Outputs:
%       massAtt     - mass attenuation coeffient at each photon energy ( cm^2 / g )
%       massAttEff  - effective mass attenuation coeffient ( cm^2 / g )
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
% 2013 - 2015

density = materialDensity( material );


if nargin == 2
    [mu, muEff] = materialAttenuation( energy, material );
elseif nargin == 3
    [mu, muEff] = materialAttenuation( energy, material, photonsPerEnergy );
else
    [mu, muEff] = materialAttenuation( energy, material, photonsPerEnergy, energyEff );
end


massAtt = mu / density;
massAttEff = muEff / density;

end

