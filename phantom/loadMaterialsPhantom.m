function [phan, map] = loadMaterialsPhantom( p )
% Load materials and density defined Phantom
%
% Meng Wu at Stanford University
% 2014 - 2015

fprintf(['Loading material phantom "' p.Phantom.materialsFileName '" with material mapping "' p.Phantom.materialMappingName '"... ']);

%%  read phantom matreial data file
[phantomMaterials, metaMaterials] = readMetaImage( [ p.Phantom.materialsFileName '.mhd' ] );

phan.size           = metaMaterials.DimSize; % number of pixels in ground truth materials image
phan.spacing        = metaMaterials.ElementSpacing; % input spacing (mm)
phan.offset         = p.Reconstruction.offset;
phan.reconOffset    = p.Reconstruction.offset;
phan.reconSize      = p.Reconstruction.size;
phan.reconSpacing   = p.Reconstruction.spacing;

fprintf('done.\n');

%% Load material mapping infomations

fprintf(['\t and applying material mapping "' p.Phantom.materialMappingName '"... ']);

[ materialIndexes, materialNames, mapTissueTemp, mapBoneTemp] = loadPantomMaterialInfo( p, phantomMaterials(:,:,ceil(end/2)) );

%% Change to 2D case if necessary

if length(phan.reconSize) == 2
    phan.size = phan.size(1:2);
    phan.spacing = phan.spacing(1:2);
    phan.offset = phan.offset(1:2);
    phantomMaterials = phantomMaterials(:,:,ceil(end/2));
end

%% Load phantom materials

phan.materialIndexes        = materialIndexes;
phan.materialNames          = materialNames;
phan.materialMappingName    = p.Phantom.materialMappingName;
phan.phantomMaterials       = phantomMaterials;
phan.materialsDir           = p.Paths.materialsDir;
phan.useDensity             = 0;


%% Scale to tissue and bone map reconstruction coordiantes

y = ( -(phan.size(1)-1)/2 : (phan.size(1)-1)/2  + phan.offset(1) ) * phan.spacing(1);
x = ( -(phan.size(2)-1)/2 : (phan.size(2)-1)/2  + phan.offset(2) )  * phan.spacing(2);

[xx, yy] = meshgrid(x, y);

iy = ( -(phan.reconSize(1)-1)/2 : (phan.reconSize(1)-1)/2  + phan.reconOffset(1) ) * phan.reconSpacing(1) ;
ix = ( -(phan.reconSize(2)-1)/2 : (phan.reconSize(2)-1)/2  + phan.reconOffset(2) ) * phan.reconSpacing(2) ;

[ixx, iyy] = meshgrid(ix, iy);

map.mapTissue   = interp2( xx, yy, mapTissueTemp, ixx, iyy );
map.mapBone     = interp2( xx, yy, mapBoneTemp, ixx, iyy );
map.mapTissue   = ( map.mapTissue > 0.9 );
map.mapBone     = ( map.mapBone > 0.9 );

if size(map.mapTissue,1) >= 256
    se1 = strel('disk',7);
else
    se1 = strel('disk',5);
end

map.mapTissue = imerode(map.mapTissue,se1);

map.windowHu    = p.Visualization.windowHu;
map.windowAtt   = p.Visualization.windowAtt;

fprintf('done.\n');

end