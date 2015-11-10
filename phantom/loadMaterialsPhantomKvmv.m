function [phan, map, roi] = loadMaterialsPhantomKvmv(p )
% Load Materials Phantom
%
% Copyright (c) 2010-2012 by Andreas Keil, Stanford University.
% Modified by Meng Wu at 2012.9

fprintf(['Loading material phantom "' p.Phantom.materialsFileName '" with material mapping "' p.Phantom.materialMappingName '"... ']);

% read phantom data file
[phantomMaterials, metaMaterials] = readMetaImage( [ p.Phantom.materialsFileName '.mhd' ] );

phan.size           = metaMaterials.DimSize; % number of pixels in ground truth materials image
phan.spacing        = metaMaterials.ElementSpacing; % input spacing (mm)
phan.offset         = p.Reconstruction.offset;
phan.reconOffset    = p.Reconstruction.offset;
phan.reconSize      = p.Reconstruction.size;
phan.reconSpacing   = p.Reconstruction.spacing;

% 2D case
if length(phan.reconSize) == 2
    phan.size = phan.size(1:2);
    phan.spacing = phan.spacing(1:2);
    phan.offset = phan.offset(1:2);
    phantomMaterials = phantomMaterials(:,:,round(end/2));
end


% origin coordinate of the reconstruction (world coordinate of the first voxel in mm)
originRecon = -(p.Reconstruction.size-1)/2 .* p.Reconstruction.spacing;

% define materials whose attenuation profiles should be read from disk and mapped to material phantom indexes
fprintf(['and applying material mapping "' p.Phantom.materialMappingName '"... ']);

%% Load phantom ROI
if strcmpi(p.Phantom.materialsFileName, 'headMaterials2Patterns2Fillings-1120x1200')
    
    % ROI1 around the soft tissue pattern between the fillings
    roi1LT = floor(([-10 -48] - originRecon(1:2)) ./ phan.reconSpacing(1:2)) + 1; % left and top coordinate, converted from mm to voxels (1-based in Matlab!)
    roi1RB = ceil(([10 -32] - originRecon(1:2)) ./ phan.reconSpacing(1:2)) + 1; % right and bottom coordinate, converted from mm to voxels (1-based in Matlab!)
    ssimRoi1 = [roi1LT(1) roi1RB(1) roi1LT(2) roi1RB(2)];
    
    % ROI2 around the soft tissue pattern that is mostly unaffected by the metal
    roi2LT = floor(([-10 0] - originRecon(1:2)) ./ phan.reconSpacing(1:2)) + 1; % left and top coordinate, converted from mm to voxels (1-based in Matlab!)
    roi2RB = ceil(([10 16] - originRecon(1:2)) ./ phan.reconSpacing(1:2)) + 1; % right and bottom coordinate, converted from mm to voxels (1-based in Matlab!)
    ssimRoi2 = [roi2LT(1) roi2RB(1) roi2LT(2) roi2RB(2)];
    
    
elseif strcmpi(p.Phantom.materialsFileName, 'headRealistic') ||  strcmpi(p.Phantom.materialsFileName, 'headRealistic-3d')  ||  strcmpi(p.Phantom.materialsFileName, 'headRealistic2-3d')
    
    % tissue pattern size
    patternSize = [15 14];
    
    % ROI1 around the soft tissue pattern between the fillings
    % 	roi1Center = ([98, 55] + [51, 38.5]) / 2;
    roi1Center = [0 -46.5];
    roi1LT = roi1Center - 0.5*patternSize;
    roi1RB = roi1Center + 0.5*patternSize;
    roi1LT = round((roi1LT - originRecon(1:2)) ./ phan.reconSpacing(1:2)); % left and top coordinate, converted from mm to voxels (1-based in Matlab!)
    roi1RB = round((roi1RB - originRecon(1:2)) ./ phan.reconSpacing(1:2)); % right and bottom coordinate, converted from mm to voxels (1-based in Matlab!)
    ssimRoi1 = [roi1LT(1) roi1RB(1) roi1LT(2) roi1RB(2)];
    
    % ROI2 around the soft tissue pattern that is mostly unaffected by the metal
    % 	roi2Center = [69, 126];
    roi2Center = [0 36.5];
    roi2LT = roi2Center - 0.5*patternSize;
    roi2RB = roi2Center + 0.5*patternSize;
    roi2LT = round((roi2LT - originRecon(1:2)) ./ phan.reconSpacing(1:2)); % left and top coordinate, converted from mm to voxels (1-based in Matlab!)
    roi2RB = round((roi2RB - originRecon(1:2)) ./ phan.reconSpacing(1:2)); % right and bottom coordinate, converted from mm to voxels (1-based in Matlab!)
    ssimRoi2 = [roi2LT(1), roi2RB(1), roi2LT(2), roi2RB(2)];
    
elseif strcmpi(p.Phantom.materialsFileName, 'HipImplantsRealistic') ||  strcmpi(p.Phantom.materialsFileName, 'HipImplantsRealistic-3d')
    
    % tissue pattern size
    patternSize = [24 24];
    
    % ROI1 around the soft tissue pattern between the fillings
    % 	roi1Center = ([98, 55] + [51, 38.5]) / 2;
    roi1Center = [5 -15];
    roi1LT = roi1Center - 0.5*patternSize;
    roi1RB = roi1Center + 0.5*patternSize;
    roi1LT = round((roi1LT - originRecon(1:2)) ./ phan.reconSpacing(1:2)); % left and top coordinate, converted from mm to voxels (1-based in Matlab!)
    roi1RB = round((roi1RB - originRecon(1:2)) ./ phan.reconSpacing(1:2)); % right and bottom coordinate, converted from mm to voxels (1-based in Matlab!)
    ssimRoi1 = [roi1LT(1) roi1RB(1) roi1LT(2) roi1RB(2)];
    
    % ROI2 around the soft tissue pattern that is mostly unaffected by the metal
    % 	roi2Center = [69, 126];
    roi2Center = [5 45];
    roi2LT = roi2Center - 0.5*patternSize;
    roi2RB = roi2Center + 0.5*patternSize;
    roi2LT = round((roi2LT - originRecon(1:2)) ./ phan.reconSpacing(1:2)); % left and top coordinate, converted from mm to voxels (1-based in Matlab!)
    roi2RB = round((roi2RB - originRecon(1:2)) ./ phan.reconSpacing(1:2)); % right and bottom coordinate, converted from mm to voxels (1-based in Matlab!)
    ssimRoi2 = [roi2LT(1), roi2RB(1), roi2LT(2), roi2RB(2)];
    
else
    % tissue pattern size
    patternSize = [10 10];
    
    % ROI1 around the soft tissue pattern between the fillings
    % 	roi1Center = ([98, 55] + [51, 38.5]) / 2;
    roi1Center = [0 0];
    roi1LT = roi1Center - 0.5*patternSize;
    roi1RB = roi1Center + 0.5*patternSize;
    roi1LT = round((roi1LT - originRecon(1:2)) ./ phan.reconSpacing(1:2)); % left and top coordinate, converted from mm to voxels (1-based in Matlab!)
    roi1RB = round((roi1RB - originRecon(1:2)) ./ phan.reconSpacing(1:2)); % right and bottom coordinate, converted from mm to voxels (1-based in Matlab!)
    ssimRoi1 = [roi1LT(1) roi1RB(1) roi1LT(2) roi1RB(2)];
    
    % ROI2 around the soft tissue pattern that is mostly unaffected by the metal
    % 	roi2Center = [69, 126];
    roi2Center = [0 0];
    roi2LT = roi2Center - 0.5*patternSize;
    roi2RB = roi2Center + 0.5*patternSize;
    roi2LT = round((roi2LT - originRecon(1:2)) ./ phan.reconSpacing(1:2)); % left and top coordinate, converted from mm to voxels (1-based in Matlab!)
    roi2RB = round((roi2RB - originRecon(1:2)) ./ phan.reconSpacing(1:2)); % right and bottom coordinate, converted from mm to voxels (1-based in Matlab!)
    ssimRoi2 = [roi2LT(1), roi2RB(1), roi2LT(2), roi2RB(2)];
    
end


%% Load phantom materials

mapTissueTemp = zeros( size( phantomMaterials), 'single' );
mapBoneTemp = zeros( size( phantomMaterials), 'single' );

if strcmpi(p.Phantom.materialMappingName, 'v3_Soft')
    
    materialIndexes = {0, 1, 2};
    materialNames = { ...
        'Air_Dry_near_sea_level', ...
        'Tissue_Soft_ICRU-44', ...
        'Tissue_Soft_ICRU-44', ...
        };
    
    mapTissueTemp( phantomMaterials == 1 |  phantomMaterials == 2 ) = 1;
    mapBoneTemp( phantomMaterials == 2 ) = 1;
    
elseif strcmpi(p.Phantom.materialMappingName, 'v5_bhc') % Andreas' more realistic phantom definition that was obtained by segmenting a real CT scan
    
    materialIndexes = {0, 1, 2, 3, 4, 5};
    materialNames = { ...
        'Air_Dry_near_sea_level', ...
        'Tissue_Soft_ICRU-44', ...
        'B-100_Bone-Equivalent_Plastic', ...
        'Aluminum'...
        'Bone_Cortical_ICRU-44', ...
        'Ca', ...
        };
    
    mapTissueTemp( phantomMaterials == 1 ) = 1;
    mapBoneTemp( phantomMaterials == 2 ) = 1;    
    
elseif strcmpi(p.Phantom.materialMappingName, 'v3_SoftBone') % Andreas' more realistic phantom definition that was obtained by segmenting a real CT scan
    
    materialIndexes = {0, 1, 2};
    materialNames = { ...
        'Air_Dry_near_sea_level', ...
        'Tissue_Soft_ICRU-44', ...
        'Bone_Cortical_ICRU-44', ...
        };
    
    mapTissueTemp( phantomMaterials == 1 ) = 1;
    mapBoneTemp( phantomMaterials == 2 ) = 1;
    
elseif strcmpi(p.Phantom.materialMappingName, 'v6_AmalgamGold')
    
    materialIndexes = {0, 1, 2, 3, 4, 5};
    materialNames = { ...
        'Air_Dry_near_sea_level', ...
        'Adipose_Tissue_ICRU-44', ...
        'Tissue_Soft_ICRU-44', ...
        'Bone_Cortical_ICRU-44', ...
        'AK-Amalgam', ... % Amalagam filling
        'Gold' ... % gold filling
        };
    
    mapTissueTemp( phantomMaterials == 1 | phantomMaterials == 2) = 1;
    mapBoneTemp( phantomMaterials == 3 ) = 1;
    
elseif strcmpi(p.Phantom.materialMappingName, 'v6_Titanium')
    
    materialIndexes = {0, 1, 2, 3, 4, 5};
    materialNames = { ...
        'Air_Dry_near_sea_level', ...
        'Adipose_Tissue_ICRU-44', ...
        'Tissue_Soft_ICRU-44', ...
        'B-100_Bone-Equivalent_Plastic', ...
        'Bone_Cortical_ICRU-44', ...
        'Titanium', ... % large filling filling
        };
    
    mapTissueTemp( phantomMaterials == 2 ) = 1;
    mapBoneTemp( phantomMaterials == 3 | phantomMaterials == 4 ) = 1;
    
    
elseif strcmpi(p.Phantom.materialMappingName, 'v6_Bones')
    
    materialIndexes = {0, 1, 2, 3, 4, 5};
    materialNames = { ...
        'Air_Dry_near_sea_level', ...
        'Adipose_Tissue_ICRU-44', ...
        'Tissue_Soft_ICRU-44', ...
        'B-100_Bone-Equivalent_Plastic', ...
        'Bone_Cortical_ICRU-44', ...
        'Bone_Cortical_ICRU-44', ...
        };
    
    mapTissueTemp( phantomMaterials == 2 ) = 1;
    mapBoneTemp( phantomMaterials == 3 | phantomMaterials == 4 ) = 1;
    
elseif strcmpi(p.Phantom.materialMappingName, 'v6_Nofillings')
    materialIndexes = {0, 1, 2, 3, 4, 5};
    
    materialNames = { ...
        'Air_Dry_near_sea_level', ...
        'Adipose_Tissue_ICRU-44', ...
        'Tissue_Soft_ICRU-44', ...
        'Bone_Cortical_ICRU-44', ...
        'Bone_Cortical_ICRU-44', ... % replace Amalagam filling with tooth
        'Bone_Cortical_ICRU-44' ... % replace gold filling with tooth
        };
    
    mapTissueTemp( phantomMaterials == 1 | phantomMaterials == 2) = 1;
    mapBoneTemp( phantomMaterials == 3 | phantomMaterials == 4 | phantomMaterials == 5 ) = 1;
    
elseif strcmpi(p.Phantom.materialMappingName, 'v7_AmalgamGold')
    
    materialIndexes = {0, 1, 2, 3, 4, 5, 6};
    materialNames = { ...
        'Air_Dry_near_sea_level', ...
        'Adipose_Tissue_ICRU-44', ...
        'Tissue_Soft_ICRU-44', ...
        'B-100_Bone-Equivalent_Plastic', ...
        'Bone_Cortical_ICRU-44', ...
        'AK-Amalgam', ... % large filling filling
        'Gold' ... % small filling
        };
    
    mapTissueTemp(  phantomMaterials == 2) = 1;
    mapBoneTemp( phantomMaterials == 3 | phantomMaterials == 4 ) = 1;
    
 elseif strcmpi(p.Phantom.materialMappingName, 'v7_Amalgam')
    
    materialIndexes = {0, 1, 2, 3, 4, 5, 6};
    materialNames = { ...
        'Air_Dry_near_sea_level', ...
        'Adipose_Tissue_ICRU-44', ...
        'Tissue_Soft_ICRU-44', ...
        'B-100_Bone-Equivalent_Plastic', ...
        'Bone_Cortical_ICRU-44', ...
        'AK-Amalgam', ... % large filling filling
        'AK-Amalgam' ... % small filling
        };
    
    mapTissueTemp(  phantomMaterials == 2) = 1;
    mapBoneTemp( phantomMaterials == 3 | phantomMaterials == 4 ) = 1;   
    
elseif strcmpi(p.Phantom.materialMappingName, 'v7_Gold')
    materialIndexes = {0, 1, 2, 3, 4, 5, 6};
    materialNames = { ...
        'Air_Dry_near_sea_level', ...
        'Adipose_Tissue_ICRU-44', ...
        'Tissue_Soft_ICRU-44', ...
        'B-100_Bone-Equivalent_Plastic', ...
        'Bone_Cortical_ICRU-44', ...
        'Gold', ... % large filling filling
        'Gold' ... % small filling
        };
    
    mapTissueTemp(  phantomMaterials == 2) = 1;
    mapBoneTemp( phantomMaterials == 3 | phantomMaterials == 4 ) = 1;
    
    
elseif strcmpi(p.Phantom.materialMappingName, 'v7_Nofillings')
    
    materialIndexes = {0, 1, 2, 3, 4, 5, 6};
    materialNames = { ...
        'Air_Dry_near_sea_level', ...
        'Adipose_Tissue_ICRU-44', ...
        'Tissue_Soft_ICRU-44', ...
        'B-100_Bone-Equivalent_Plastic', ...
        'Bone_Cortical_ICRU-44', ...
        'Bone_Cortical_ICRU-44', ...
        'Bone_Cortical_ICRU-44'
        };
    
    mapTissueTemp(  phantomMaterials == 2) = 1;
    mapBoneTemp( phantomMaterials == 3 | phantomMaterials == 4 ) = 1;
    
elseif strcmpi(p.Phantom.materialMappingName, 'p5_Iodine') % Andreas' more realistic phantom definition that was obtained by segmenting a real CT scan
    
    materialIndexes = {0, 1, 2, 3, 4};
    materialNames = { ...
        'air', ...
        'water', ...
        'brain', ...
        'bone_cortical', ...
        'I', ...
        };
    
    mapTissueTemp( phantomMaterials == 1 ) = 1;
    mapBoneTemp( phantomMaterials == 2 ) = 1;
    
else
    
    error('Unknown phantom type selected!');
    
end


phan.materialIndexes        = materialIndexes;
phan.materialNames          = materialNames;
phan.materialMappingName    = p.Phantom.materialMappingName;
phan.phantomMaterials       = phantomMaterials;
phan.materialsDir           = p.Paths.materialsDir;
phan.useDensity             = 0;

%% Scale to reconstruction

y = [ -(phan.size(1)-1)/2: (phan.size(1)-1)/2] * phan.spacing(1) + phan.offset(1);
x = [ -(phan.size(2)-1)/2: (phan.size(2)-1)/2] * phan.spacing(2) + phan.offset(2);

[xx, yy] = meshgrid(x, y);

iy = [ -(phan.reconSize(1)-1)/2: (phan.reconSize(1)-1)/2] * phan.reconSpacing(1) + phan.reconOffset(1);
ix = [ -(phan.reconSize(2)-1)/2: (phan.reconSize(2)-1)/2] * phan.reconSpacing(2) + phan.reconOffset(2);

[ixx, iyy] = meshgrid(ix, iy);

if length( phan.size ) == 2
    
    map.mapTissue = interp2( xx, yy, mapTissueTemp, ixx, iyy );
    map.mapBone = interp2( xx, yy, mapBoneTemp, ixx, iyy );
    
else
    map.mapTissue = interp2( xx, yy, squeeze( mapTissueTemp( :,:,ceil(phan.size(3)/2))), ixx, iyy );
    map.mapBone = interp2( xx, yy, squeeze( mapBoneTemp( :,:,ceil(phan.size(3)/2))), ixx, iyy );
end

map.mapTissue   = ( map.mapTissue > 0.9 );
map.mapBone     = ( map.mapBone > 0.9 );

if size(map.mapTissue,1) > 256
    se1 = strel('disk',9);
else
    se1 = strel('disk',5);
end

map.mapTissue = imerode(map.mapTissue,se1);

roi.size        = phan.reconSize;
roi.windowHu    = p.Visualization.windowHu;
%roi.roi1Center  = roi1Center;
roi.ssimRoi1    = ssimRoi1;
%roi.roi2Center  = roi2Center;
roi.ssimRoi2    = ssimRoi2;

fprintf('done.\n');
fprintf('\n');


end