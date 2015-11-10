function img = reconHelicalWeightedFDK( sino, geom, H, noViewsPerSlice, tiltedFilter )
% function img = reconHelicalWeightedFDK( sino, geom, window, segmentLength, tiltedFilter, crop )
% FDK helical CT recontruction not efficient in weighting and filtering, I
% guess it more robust in noisy data than the other version
% input:
%       sino    - log sinogram
%       geom    - geometry parameters
%       window  - window function ( default 'hamming', other options: 'ram-lak', 'shepp-logan', 'cosine'  )
%       viewUpsamplingRate (default 1)
%       crop    - frequency crop ratio (default 1)
%       mode    - length of segement to reconstruct a slice ( default 0: pi + fan scan)
% output:
%       img     - reconstructed attenuation image (1/cm)
%
% The first FDK algorith in Kudo et al, "Exact and appoximate algorithm for helical
% cone-beam CT," PMB, 2004
%
% Meng Wu, Stanford University, 2014

if nargin < 5,  tiltedFilter = false;   end

nz          = geom.reconSize(3);
oz          = geom.reconOffset(3);
dz          = geom.reconSpacing(3);
noViews     = geom.noViews;
couchZ      = geom.couchZ;

dbeta = abs(geom.betas(end) - geom.betas(1)) / (noViews-1);

% compute redundant weighted for one slice reconstruction
redundantWeights = computeHelicalWeights2D( geom, noViewsPerSlice );

geomSlice = geom;
geomSlice.reconSize(3)  = 1;
geomSlice.helicalScan   = 0;

% cosine weighting
sino = cosineWeighting( sino, geom );

img = zeros( geom.reconSize, 'single' );

fprintf( '\tSlice: ');
for iz = 1: nz
    
    if mod(iz, 10 ) == 0 || iz == 1
        fprintf( '%i/(%i), ', iz, nz );
    end
    
    % slice position to reconstruct
    z = ( iz - (nz+1)/2 + oz ) * dz ;
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
    
    sinoSlice = filterHelicalSinogram( sino(:,:,iview_start:iview_stop), H, redundantWeights, geom, tiltedFilter );
    
    % finally, single slice FBP recon
    img(:,:,nz-iz+1) = backProjectMex(sinoSlice, geomSlice,  1, 0, 'back,hpd' ) * dbeta * 10;
end


end



