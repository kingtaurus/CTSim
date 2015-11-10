function [ztable, density ] = findCompositionMaterial( material )

if nargin == 0
    ztable = {'adipose', 'blood', 'bone_compact', 'bone_cortical', 'brain', 'lung', ...
        'muscle_skeletal', 'muscle_striated', 'skin', 'soft_tissue', 'water', 'air' ...
        'CWO', 'Acrylic', 'PMP', 'Delrin', 'Teflon', 'Polystryrene', 'Bone20', 'Bone50' ...,
        'LDPE', 'Quartz' };
    return;
end

compositionData = BodyMaterialCompositionFunc;

switch material
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
    case 'air'
        ztable = compositionData.air(2:end,:);
        density = compositionData.air(1,2);
    case 'CWO'
        ztable = compositionData.cwo(2:end,:);
        density = compositionData.cwo(1,2);
    case 'Acrylic'
        ztable = compositionData.acrylic(2:end,:);
        density = compositionData.acrylic(1,2);
    case 'PMP'
        ztable = compositionData.pmp(2:end,:);
        density = compositionData.pmp(1,2);
    case 'Delrin'
        ztable = compositionData.delrin(2:end,:);
        density = compositionData.delrin(1,2);
    case 'Teflon'
        ztable = compositionData.teflon(2:end,:);
        density = compositionData.teflon(1,2);
    case 'Polystryrene'
        ztable = compositionData.polystryrene(2:end,:);
        density = compositionData.polystryrene(1,2);
    case 'Bone20'
        ztable = compositionData.bone20(2:end,:);
        density = compositionData.bone20(1,2);
    case 'Bone50'
        ztable = compositionData.bone50(2:end,:);
        density = compositionData.bone50(1,2);
    case 'LDPE'
        ztable = compositionData.ldpe(2:end,:);
        density = compositionData.ldpe(1,2);
    case 'Quartz'
        ztable = compositionData.quartz(2:end,:);
        density = compositionData.quartz(1,2);
    otherwise
        ztable = material;
        [ ~ , ~, element_data] = XrayData;
        z = ChemicalSymbols2AtomicNumbers({material});
        density = element_data{z}.Density;
end

end