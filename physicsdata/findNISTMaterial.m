function inNIST = findNISTMaterial( material )

list = {'Air_Dry_near_sea_level', 'Adipose_Tissue_ICRU-44', 'Water_Liquid', 'Lung_Tissue_ICRU-44', ...
    'Tissue_Soft_ICRU-44', 'Muscle_Skeletal_ICRU-44', 'B-100_Bone-Equivalent_Plastic', ...
    'Bone_Cortical_ICRU-44', 'Aluminum', 'Titanium', 'Iron', 'Copper', 'AK-Amalgam', 'Gold'};

if nargin == 0
    inNIST = list;
end

x = strcmp( material, list );

if  sum( x ) > 0
    inNIST = true;
else
    inNIST = false;
end
   
end