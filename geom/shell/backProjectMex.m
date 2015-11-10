function img = backProjectMex( sino, geom, M, i, projType, map )
% Backward projection using distance driven method with ordered subset in 2D
%       and 3d. Exact traspose to the forwardProjectDistanceDrivenMex
% input:
%       sino - projection result aka sinogram
%       geom - system geometry
%       M    - number of subsets
%       i    - subset number
% output:
%       img - image to project
%
% Meng Wu
% 2013.4

if nargin < 3
    M = 1;
    i = 0;
end

if nargin < 5
    projType = 'back,dd';
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
    
    img = fbct_geom_mex( projType, int32(geom.reconSize), single(geom.reconSpacing), ...
        single(geom.reconOffset), int32(geom.detSize), single(geom.detSpacing), ...
        single(geom.detOffset), single(geom.SAD), single(geom.ADD), single(geom.SDD), ...
        int32(noViews), double(betas),  int32(geom.flatPanel), single(sino) );
    
elseif isfield( geom, 'PMat') %carm CT
    
    img = cbct_geom_pmat_mex( projType, int32(geom.reconSize), single(geom.reconSpacing), ...
        single(geom.reconOffset), int32(geom.detSize), single( geom.detSpacing), int32(noViews), single(geom.PMat), ...
        single(geom.cameraPositions), single(geom.SDDPerProjection),  single( geom.detOffsetsPerProjection), single(sino), logical(geom.map) );
    
elseif length(geom.reconSize) == 3 && length(geom.detSize) == 2
    
    if ( strcmp( projType, 'back,pd') || strcmp( projType, 'back,hpd') ) && USE_CUDA_DEVICE
        
        img = cbct_back_pd_cuda( projType, int32(geom.reconSize), single(geom.reconSpacing), ...
            single(geom.reconOffset), int32(geom.detSize), single(geom.detSpacing), ...
            single(geom.detOffset), single(geom.SAD), single(geom.ADD), single(geom.SDD), ...
            int32(noViews), double(betas), single(couchZ), int32(geom.flatPanel), single(sino), logical(geom.map) );
        
    elseif  strcmp( projType, 'back,dd') && USE_CUDA_DEVICE
        
        img = cbct_dd_cuda( projType, int32(geom.reconSize), single(geom.reconSpacing), ...
            single(geom.reconOffset), int32(geom.detSize), single(geom.detSpacing), ...
            single(geom.detOffset), single(geom.SAD), single(geom.ADD), single(geom.SDD), ...
            int32(noViews), double(betas), single(couchZ), int32(geom.flatPanel), single(sino), logical(geom.map) );
        
    elseif  strcmp( projType, 'back,tf') && USE_CUDA_DEVICE
        
        img = cbct_sf_cuda( projType, int32(geom.reconSize), single(geom.reconSpacing), ...
            single(geom.reconOffset), int32(geom.detSize), single(geom.detSpacing), ...
            single(geom.detOffset), single(geom.SAD), single(geom.ADD), single(geom.SDD), ...
            int32(noViews), double(betas), single(couchZ), int32(geom.flatPanel), single(sino), logical(geom.map) );
        
    else
        
        img = cbct_geom_mex( projType, int32(geom.reconSize), single(geom.reconSpacing), ...
            single(geom.reconOffset), int32(geom.detSize), single(geom.detSpacing), ...
            single(geom.detOffset), single(geom.SAD), single(geom.ADD), single(geom.SDD), ...
            int32(noViews), double(betas), single(couchZ), int32(geom.flatPanel), single(sino), logical(geom.map) );
    end
    
else
    Error('Wrong geometry dimension!\n');
end

if M ~= 1
    img = M *img;
end

end