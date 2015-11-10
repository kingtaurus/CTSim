function img = reconHelicalFBP( sino, geom, window, noViewsPerSlice, tiltedFilter, crop )
% function img = reconHelicalFBP( sino, geom, window, segmentLength, tiltedFilter, crop )
% Filtered back-projection of helical CT recontruction using linear
% interpolation method to get sinogram for each slice.
%
% Very old technique, only work for the narrow detector and small pitch
%
% input:
%       sino    - log sinogram
%       geom    - geometry parameters
%       window  - window function ( default 'hamming', other options: 'ram-lak', 'shepp-logan', 'cosine'  )
%       segmentLength
%       tiltedFilter
%       crop    - frequency crop ratio (default 1)
% output:
%       img     - reconstructed attenuation image (1/cm)
%
% Meng Wu, Stanford University, 2014

if nargin < 5,  tiltedFilter = false;   end
if nargin < 6,  crop = 1;   end

nz = geom.reconSize(3);
oz = geom.reconOffset(3);
dz = geom.reconSpacing(3);
couchZ      = geom.couchZ;

u = ( ( -(geom.detSize(1)-1)/2:(geom.detSize(1)-1)/2) + geom.detOffset(1) ) * geom.detSpacing(1);
v = ( ( -(geom.detSize(2)-1)/2:(geom.detSize(2)-1)/2) + geom.detOffset(2) ) * geom.detSpacing(2);
[uu, vv ] = meshgrid( u, v);

% weighting factor and ramp filter
if geom.flatPanel
    gammas   = atan( - u / geom.SDD );
else
    gammas   = - u / geom.SDD ;
end

dbeta = abs(geom.betas(end) - geom.betas(1)) / (geom.noViews-1);
vt   = gammas / dbeta * geom.couchSpacing;

if geom.detSize(2) > 8
    fprintf('\tWarning: FBP method is not suitable for large cone angle. \n');
end

% parameter for single slice FBP recon
geomSlice = geom;
geomSlice.detSize       = geomSlice.detSize(1);
geomSlice.detSpacing    = geomSlice.detSpacing(1);
geomSlice.detOffset     = geomSlice.detOffset(1);
geomSlice.reconSize     = geomSlice.reconSize(1:2);
geomSlice.reconSpacing  = geomSlice.reconSpacing(1:2);
geomSlice.reconOffset   = geomSlice.reconOffset(1:2);
geomSlice.helicalScan   = 0;

% cosine weighting
sino = cosineWeighting( sino, geom );
img = zeros( geom.reconSize, 'single' );

fprintf('\tSlice: ');
for iz = 1: nz
    
    if mod(iz, 10 ) == 0
        fprintf( '%i/(%i), ', iz, nz );
    end
    
    % slice position to reconstruct
    z = ( iz - (nz+1)/2 + oz ) * dz ;
    [ ~, iview_middle ] = min(abs(couchZ-z));
    
    % center source postion
    if couchZ(iview_middle) - z > 0
        iview_middle = iview_middle - 1;
    end
    
    % determine the begin and end veiw
    iview_start = iview_middle - noViewsPerSlice/2 + 1;
    iview_stop  = iview_middle + noViewsPerSlice/2;
    
    iview_start = max( iview_start, 1 );
    iview_stop  = min( iview_stop, geom.noViews );
    
    % not enough projection to do a short scan so just scape this slice
    if iview_stop - iview_start < geom.noViewsTurn /2
        continue;
    end
    
    geomSlice.noViews = iview_stop - iview_start + 1;
    geomSlice.betas = geom.betas( iview_start : iview_stop );
    
    
    if geomSlice.noViews == geom.noViewsTurn
        geomSlice.shortScan = 0;
    else
        geomSlice.shortScan = 1;
    end
    
    sinoSlice = zeros( [geom.detSize(1) geomSlice.noViews], 'single' );
    
    %linear interpolate sinogram for a slice
    for iview = 1:geomSlice.noViews
        
        % z position on detecor
        v_slice = couchZ( iview + iview_start - 1 ) - z;
        
        if tiltedFilter
            slice =  interp2( uu, vv, squeeze( sino(:,:,iview + iview_start - 1) ),  u, v_slice - vt  ) ;
        else
            slice =  interp2( uu, vv, squeeze( sino(:,:,iview + iview_start - 1) ),  u, v_slice ) ;
        end
        
        slice( isnan(slice) ) = 0;
        sinoSlice(:, iview) = slice;
        
    end
    
    % finally, single slice FBP recon
    img(:,:,nz-iz+1) =  reconFBP2d( sinoSlice, geomSlice, window, crop );
    
end

if geomSlice.shortScan
    img = img * 2;
end

