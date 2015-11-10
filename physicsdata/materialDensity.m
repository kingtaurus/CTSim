function density = materialDensity( material )
% Return density of materails at given x-ray energy ( g / cm^2)
% Inputs: 
%       material - material name, see below for defined material names
% Outputs:
%       density - ( g / cm^2)
%
% NIST materials:
%       {'Air_Dry_near_sea_level', 'Adipose_Tissue_ICRU-44', 'Water_Liquid', ...
%       'Tissue_Soft_ICRU-44', 'Muscle_Skeletal_ICRU-44', 'B-100_Bone-Equivalent_Plastic', ...
%       'Bone_Cortical_ICRU-44', 'Aluminum', 'Titanium', 'Iron', 'Copper', 'AK-Amalgam', 'Gold'}
% Self defined materials:
%       {'adipose', 'blood', 'bone_compact', 'bone_cortical', 'brain', 'lung', ...
%        'muscle_skeletal', 'muscle_striated', 'skin', 'soft_tissue', 'water', 'air' ...
%        'CWO', 'Acrylic', 'PMP', 'Delrin', 'Teflon', 'Polystryrene', 'Bone20', 'Bone50' ...,
%        'LDPE', 'Quartz' };
% Elements with alphabetical symbles
%
% Meng Wu at Stanford University
% 2013 - 2014

if findNISTMaterial( material )
    [~, ~, density] = readMaterialAttenuation(material, 'physicsdata/materials/');
else
    [~, density ] = findCompositionMaterial( material );
end


end

