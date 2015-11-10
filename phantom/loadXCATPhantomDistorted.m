function [phan, map ] = loadXCATPhantomDistorted( p )
% Load Materials Phantom
%
% Copyright (c) 2010-2012 by Andreas Keil, Stanford University.
% Modified by Meng Wu at 2012.9
% Modified by He Yang at 2015.1

fprintf(['Loading material phantom "' p.Phantom.materialsFileName '" with material mapping "' p.Phantom.materialMappingName '"... ']);
%fprintf(['Loading material phantom "' XCATPelvis1 '" with material mapping "' p.Phantom.materialMappingName '"... ']);

% read phantom data file from CTData/PhantomData/xcat
%[phantomMaterials, metaMaterials] = readMetaImage( [ p.Phantom.materialsFileName '.mhd' ] );
%[phantomDensities ] = readMetaImage( [ p.Phantom.materialsFileName '-density.mhd' ] );
[phantomMaterials, metaMaterials] = readMetaImage(  'XCATPelvis2-3d.mhd'  );
[phantomDensities ] = readMetaImage(  'XCATPelvis2-3d-density.mhd'  );
phan.size           = metaMaterials.DimSize; % number of pixels in ground truth materials image
phan.spacing        = metaMaterials.ElementSpacing; % input spacing (mm)
phan.offset         = p.Reconstruction.offset;


phan.reconSize      = p.Reconstruction.size;
phan.reconSpacing   = p.Reconstruction.spacing;
phan.reconOffset    = p.Reconstruction.offset;

% 2D case
if length(phan.reconSize) == 2
    phan.size = phan.size(1:2);
    phan.spacing = phan.spacing(1:2);
    phan.offset = phan.offset(1:2);
    phantomMaterials = phantomMaterials(:,:,round(end/2));
    phantomDensities = phantomDensities(:,:,round(end/2));
end


% define materials whose attenuation profiles should be read from disk and mapped to material phantom indexes
fprintf(['and applying material mapping "' p.Phantom.materialMappingName '"... ']);


%% Load phantom materials

mapTissueTemp = zeros( size( phantomMaterials), 'single' );
mapBoneTemp = zeros( size( phantomMaterials), 'single' );

if strcmpi(p.Phantom.materialMappingName, 'Catphan-sensitivity')
    
    materialIndexes = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14 };
    materialNames = { ...
        'water', ...
        'Acrylic', ...
        'Bone50', ...
        'Acrylic', ...
        'PMP', ...
        'air', ...
        'lung', ...
        'Delrin', ...
        'Polystryrene', ...
        'Teflon', ...
        'Bone20', ...
        'LDPE', ...
        'Acrylic', ...
        };
elseif strcmpi(p.Phantom.materialMappingName, 'Catphan-contrast')
    
    materialIndexes = { 1 };
    materialNames = { ...
        'water' ...
        };
elseif strcmpi(p.Phantom.materialMappingName, 'Catphan-geometry')
    
    materialIndexes = { 1, 100 };
    materialNames = { ...
        'water', ...
        'Iron' ...
        };
elseif strcmpi(p.Phantom.materialMappingName, 'Catphan-resolution')
    
    materialIndexes = { 1, 2 };
    materialNames = { ...
        'water', ...
        'Titanium' ...
        };    
elseif strcmpi(p.Phantom.materialMappingName, 'Catphan-resolution-x2')
    
    materialIndexes = { 1, 2 };
    materialNames = { ...
        'water', ...
        'bone_cortical' ...
        };  
    
    
elseif strcmpi(p.Phantom.materialMappingName, 'v2-water-bone')
    
    materialIndexes = { 1, 2 };
    materialNames = { ...
        'water', ...
        'bone_compact' ...
        };
elseif strcmpi(p.Phantom.materialMappingName, 'v4-XCAT-lung')
    
    materialIndexes = {1, 2, 3, 4};
    materialNames = { ...
        'Lung_Tissue_ICRU-44', ...
        'Tissue_Soft_ICRU-44', ...
        'B-100_Bone-Equivalent_Plastic', ...
        'Bone_Cortical_ICRU-44' ...
        };
    
    mapTissueTemp( phantomMaterials == 2 ) = 1;
    mapBoneTemp( phantomMaterials == 3 | phantomMaterials == 4 ) = 1;
    
elseif strcmpi(p.Phantom.materialMappingName, 'v4-XCAT-head')
    
    materialIndexes = {1, 2, 3, 4};
    materialNames = { ...
        'soft_tissue', ...
        'brain', ...
        'muscle_striated', ...
        'bone_cortical' ...
        };
    
    mapTissueTemp( phantomMaterials == 2 ) = 1;
    mapBoneTemp( phantomMaterials == 4 ) = 1;    
elseif strcmpi(p.Phantom.materialMappingName, 'v6-XCAT-liver')
    
    materialIndexes = {1, 2, 3, 4, 5, 6};
    materialNames = { ...
        'lung', ...
        'soft_tissue', ...
        'muscle_skeletal', ...
        'blood' ...
        'bone_compact', ...
        'bone_cortical' ...
        };
    
    mapTissueTemp( phantomMaterials == 4 ) = 1;
    mapBoneTemp( phantomMaterials == 5 | phantomMaterials == 6 ) = 1;
    
    

else
    
    error('Unknown phantom type selected!');
    
end


phan.materialIndexes        = materialIndexes;
phan.materialNames          = materialNames;
phan.materialMappingName    = p.Phantom.materialMappingName;
phan.phantomMaterials       = phantomMaterials;
phan.phantomDensities       = phantomDensities;
phan.materialsDir           = p.Paths.materialsDir;
phan.useDensity             = 1;

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
    se1 = strel('disk',7);
else
    se1 = strel('disk',5);
end

map.mapTissue = imerode(map.mapTissue,se1);

map.windowHu    = p.Visualization.windowHu;
map.windowAtt   = p.Visualization.windowAtt;

fprintf('done.\n');
fprintf('\n');


end