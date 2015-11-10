function sino = forwardProjectMex( img, geom, M, i,  projType, map )
% Forward projection
% input:
%       img - image to project
%       geom - system geometry
%       M    - number of subsets
%       i    - subset number
% output:
%       sinp - projection result aka sinogram
%
% Meng Wu
% 2011.12

if nargin < 3
    M = 1;
    i = 0;
end

if nargin < 5
    projType = 'proj,dd';
end

if nargin >= 6
    geom.map = map;
end

if ~isfield( geom, 'map' )
    geom.map = true( [geom.reconSize(1) geom.reconSize(2)] );
end

if M == 1
    betas = geom.betas;
    noViews = geom.noViews;
    couchZ  = geom.couchZ;
else
    betas   = geom.betas(1+mod(i,M):M:end);
    noViews = length(betas);
    couchZ  = geom.couchZ(1+mod(i,M):M:end);
end

if license('test','Distrib_Computing_Toolbox') && ispc
    g = gpuDevice;
    USE_CUDA_DEVICE = g(1).ToolkitVersion >= 5 ;
else
    USE_CUDA_DEVICE = false;
end



if length(geom.reconSize) == 2 && length(geom.detSize) == 1
    
    sino = fbct_geom_mex( projType, int32(geom.reconSize), single(geom.reconSpacing), ...
        single(geom.reconOffset), int32(geom.detSize), single(geom.detSpacing), ...
        single(geom.detOffset), single(geom.SAD), single(geom.ADD), single(geom.SDD), ...
        int32(noViews), double(betas),  int32(geom.flatPanel), single(img) );

elseif isfield( geom, 'PMat') %carm CT
    
    sino = cbct_geom_pmat_mex( 'proj,dd', int32(geom.reconSize), single(geom.reconSpacing), ...
        single(geom.reconOffset), int32(geom.detSize), single( geom.detSpacing), int32(noViews), single(geom.PMat), ...
        single(geom.cameraPositions), single(geom.SDDPerProjection),  single( geom.detOffsetsPerProjection), single(img), logical(geom.map) );
    
    
elseif  length(geom.reconSize) == 3 && length(geom.detSize) == 2
    
    if strcmp( projType, 'proj,dd') && USE_CUDA_DEVICE
        
        sino = cbct_dd_cuda( projType, int32(geom.reconSize), single(geom.reconSpacing), ...
            single(geom.reconOffset), int32(geom.detSize), single(geom.detSpacing), ...
            single(geom.detOffset), single(geom.SAD), single(geom.ADD), single(geom.SDD), ...
            int32(noViews), double(betas), single(couchZ), int32(geom.flatPanel), single(img), logical(geom.map) );
        
    elseif strcmp( projType, 'proj,tf') && USE_CUDA_DEVICE
        
        sino = cbct_sf_cuda( projType, int32(geom.reconSize), single(geom.reconSpacing), ...
            single(geom.reconOffset), int32(geom.detSize), single(geom.detSpacing), ...
            single(geom.detOffset), single(geom.SAD), single(geom.ADD), single(geom.SDD), ...
            int32(noViews), double(betas), single(couchZ), int32(geom.flatPanel), single(img), logical(geom.map) );
        
    else
        
        sino = cbct_geom_mex( projType, int32(geom.reconSize), single(geom.reconSpacing), ...
            single(geom.reconOffset), int32(geom.detSize), single(geom.detSpacing), ...
            single(geom.detOffset), single(geom.SAD), single(geom.ADD), single(geom.SDD), ...
            int32(noViews), double(betas), single(couchZ), int32(geom.flatPanel), single(img), logical(geom.map) );
        
    end
    
else
    error('Wrong geometry dimension!\n');
end



end