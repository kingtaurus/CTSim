%function [ mus ] = body_material_mu( bodyMaterial, egys, compositionData )
% Return X-ray attenuation coefficients for common bodymaterial
% input: 
%       bodyMaterial - type of common body material including
%         adipose
%         blood
%         bone_compact  
%         bone_cortical
%         brain
%         lung
%         muscle_skeletal
%         muscle_striated
%         skin
%         soft_tissue
%         water
%       egys -  vector of x-ray energy (keV) 
%       compositionData -  a struct with fields composition tables for some body materials
% ouput: 
%       mus - attenuation coefficients at the specified energies ( mm^-1 )
%
% Based on John D'Errico's (2009) code at http://aprendtech.com/wordpress/?p=45
% Meng Wu 
% 2012.1
%


function [ mus ] = body_material_mu( bodyMaterial, egys, compositionData )

if nargin < 3
    compositionData = BodyMaterialCompositionFunc;
end


switch bodyMaterial
    case 'adipose'
        ztable = compositionData.adipose(2:end,:);
        density = compositionData.adipose(1,2);
    case 'blood'
        ztable = compositionData.blood(2:end,:);
        density = compositionData.blood(1,2);
    case 'bone_compact'
        ztable = compositionData.bone_compact(2:end,:);
        density = compositionData.bone_compact(1,2);
    case 'bone_cortical'
        ztable = compositionData.bone_cortical(2:end,:);
        density = compositionData.bone_cortical(1,2);
    case 'brain'
        ztable = compositionData.brain(2:end,:);
        density = compositionData.brain(1,2);
    case 'lung'
        ztable = compositionData.lung(2:end,:);
        density = compositionData.lung(1,2);
    case 'muscle_skeletal'
        ztable = compositionData.muscle_skeletal(2:end,:);
        density = compositionData.muscle_skeletal(1,2);
    case 'muscle_striated'
        ztable = compositionData.muscle_striated(2:end,:);
        density = compositionData.muscle_striated(1,2);
    case 'skin'
        ztable = compositionData.skin(2:end,:);
        density = compositionData.skin(1,2);
    case 'soft_tissue'
        ztable = compositionData.soft_tissue(2:end,:);
        density = compositionData.soft_tissue(1,2);
    case 'water'
        ztable = compositionData.water(2:end,:);
        density = compositionData.water(1,2);
    otherwise
        fprintf('Error: bodyMaterial not recognized. \n');
        return;
end

[mus] = density * XrayMu(ztable,egys) ;


end