function img = reconHelicalImageWeightedFDK( sino, geom, window, segmentLength, tiltedFilter, crop )
% function img = reconHelicalImageWeightedFDK( sino, geom, window, segmentLength, tiltedFilter, crop )
% Filtered back-projection of helical CT recontruction
% input:
%       sino    - log sinogram
%       geom    - geometry parameters
%       window  - window function ( default 'hamming', other options: 'ram-lak', 'shepp-logan', 'cosine'  )
%       viewUpsamplingRate (default 1)
%       crop    - frequency crop ratio (default 1)
% output:
%       img     - reconstructed attenuation image (1/cm)
%
% There is a serious problem in correcting the cone angle weighting.
% Correct way of cone angle weighting see function
% reconHelicalViewWeightedFDK() 
% This method is abandent. 
%   
%
% Meng Wu, Stanford University, 2014


if nargin < 3,  window = 'hamming';   end
if nargin < 4,  segmentLength = 'short';    end
if nargin < 5,  tiltedFilter = false;   end
if nargin < 6,  crop = 1;   end

nz          = geom.reconSize(3);
oz          = geom.reconOffset(3);
dz          = geom.reconSpacing(3);

nu          = double( geom.detSize(1) );
nv          = double( geom.detSize(2) );
du          = geom.detSpacing(1);
dv          = geom.detSpacing(2);
ou          = geom.detOffset(1);
ov          = geom.detOffset(2);
SDD         = geom.SDD;
noViewsTurn = geom.noViewsTurn;
noViews     = geom.noViews;
h           = geom.couchSpacing;
couchZ      = geom.couchZ;

% array of detector positions in mm
u = ( ( -(nu-1)/2:(nu-1)/2) + ou ) * du;
v = ( ( -(nv-1)/2:(nv-1)/2) + ov ) * dv ;
[uu, vv ] = meshgrid( u, v);

if ~strcmpi(segmentLength, '2pi')  % 2 pi segments
    fprintf('Warnning: unknown segment length, this method only works for 2pi segment.\n');
end

% weighting factor and ramp filter
if geom.flatPanel
    H       = du*designFilter2(du, window, nu, crop);
    gammas   = atan( - u / SDD );
    tanAlpha = abs( vv ./ sqrt(  uu.^2 + vv.^2 + SDD^2 ) );
else
    H       = du*designEquiangularFilter2(du, SDD, window, nu, crop);
    gammas   = - u / SDD ;
    tanAlpha = abs( vv ./ SDD );
end

dbeta = abs(geom.betas(end) - geom.betas(1)) / (noViews-1);
vt   = vv - repmat( gammas , nv, 1) / dbeta * h;

% detect not sufficient slice coverage caused by the high pitch 
coneCoverageFactor = ( geom.SAD - geom.reconSize(1) * geom.reconSpacing(1) / 2 ) / geom.SDD ;
if geom.pitch > coneCoverageFactor
    fprintf( 'Warning: the pitch may be too large to coverage entire slice for current segment lentgh. \n');
end

% calculate number of projection in each segment
noViewsPerSlice =  noViewsTurn ;

k = 0.8  * pi / ( 2 * max( abs(tanAlpha(:)) ) );
% cached redundant weighted for a completed segment for reconstruction
cachedRedundantWeights = 0.5 * ones( [ geom.detSize(2), geom.detSize(1), noViewsPerSlice], 'single' );
for iview = 1 : noViewsPerSlice
    cachedRedundantWeights(:,:,iview) = cos( tanAlpha * k );
end

geomSlice = geom;
geomSlice.reconSize(3)  = 1;
geomSlice.helicalScan   = 0;

% cosine weighting
sino = cosineWeighting( sino, geom );
img = zeros( geom.reconSize, 'single' );
doneComputeScales = false;

for iz = 1: nz
    
    % slice position to reconstruct
    z = ( iz - (nz-1)/2 + oz ) * dz ;
    [ ~, iview_middle ] = min(abs(couchZ-z));
    
    % determine the begin and end veiw
    iview_start = iview_middle - noViewsPerSlice/2 + 1;
    iview_stop  = iview_middle + noViewsPerSlice/2;
    
    iview_start = max( iview_start, 1 );
    iview_stop  = min( iview_stop, geom.noViews );
    
    geomSlice.reconOffset(3)   = - z / dz ;
    geomSlice.noViews   = iview_stop - iview_start + 1;
    geomSlice.betas     = geom.betas(   iview_start : iview_stop );
    geomSlice.couchZ    = geom.couchZ(  iview_start : iview_stop );
    geomSlice.shortScan = 1;
    
    % not enough projection, don't want to calculate the new weights so
    % so just skip this slice
    if iview_stop - iview_start + 1 < noViewsPerSlice
        continue;
    end
    
    if ~ doneComputeScales
        rescaleWeights = backProjectMex(cachedRedundantWeights, geomSlice,  1, 0, 'back,spd' );
        rescaleWeights = 0.5 * noViewsTurn  ./ rescaleWeights;
        rescaleWeights( rescaleWeights > 2 )= 2;
        defaultAngle = geomSlice.betas(1);
        doneComputeScales = true;
    end
    
    sinoSlice = zeros( [nv nu geomSlice.noViews], 'single' );
    
    % filter the segement of sinogram for single slice reconstruction
    for iview = 1:geomSlice.noViews
        
        redundantWeights = cachedRedundantWeights(:,:,iview);
        proj = redundantWeights .* sino(:, :, iview + iview_start - 1 );
        
        if tiltedFilter
            % tilte the sinogram
            proj = tilteProjction(  proj, vv , vt );
            % ramp filter
            proj = filterFreq( proj, H, 2 ) ;
            % tilte the sinogram back
            proj = tilteProjction(  proj, vt , vv );
        else
            proj = filterFreq( proj, H, 2 ) ;
        end
        
        sinoSlice(:, :,  iview) = proj;
        
    end
    
    % finally, single slice FBP recon
    rotateRescaleWeights = imrotate( rescaleWeights,  180 * ( defaultAngle - geomSlice.betas(1) ) / pi, 'crop' );
    rotateRescaleWeights( rotateRescaleWeights == 0 ) = 0.5;
    
    img(:,:,nz-iz+1) = rotateRescaleWeights .* backProjectMex(sinoSlice, geomSlice,  1, 0, 'back,pd' ) * dbeta * 10;
    
end


end