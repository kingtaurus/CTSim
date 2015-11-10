function img = reconHelicalApproxWeightedFDK( sino, geom, H, noViewsPerSlice, tiltedFilter )
% function img = reconHelicalApproxWeightedFDK( sino, geom, window, segmentLength, tiltedFilter, crop )
% FDK helical CT recontruction a more efficient implementation using hilber
% transform
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
% The second FDK algorith in Kudo et al, "Exact and appoximate algorithm for helical
% cone-beam CT," PMB, 2004
%
% Meng Wu, Stanford University, 2014

if nargin < 5,  tiltedFilter = false;   end

if ~geom.flatPanel
    fprintf('\tWarning: the approximated method seems only work for the flat detector. \n');
end

nz          = geom.reconSize(3);
oz          = geom.reconOffset(3);
dz          = geom.reconSpacing(3);
SDD         = geom.SDD;
noViews     = geom.noViews;
couchZ      = geom.couchZ;

% array of detector positions in mm
u = ( ( -(geom.detSize(1)-1)/2:(geom.detSize(1)-1)/2) + geom.detOffset(1) ) * geom.detSpacing(1);
v = ( ( -(geom.detSize(2)-1)/2:(geom.detSize(2)-1)/2) + geom.detOffset(2) ) * geom.detSpacing(2);
[~, vv ] = meshgrid( u, v);

% weighting factor and ramp filter
if geom.flatPanel
    gammas   = atan( - u / SDD );
else
    gammas   = - u / SDD ;
end

dbeta = abs(geom.betas(end) - geom.betas(1)) / (noViews-1);
vt    = vv - repmat( gammas , geom.detSize(2), 1) / dbeta * geom.couchSpacing ;
rampFilterGain  = sum( real( fft(H)));

geomSlice = geom;
geomSlice.reconSize(3)  = 1;
geomSlice.helicalScan   = 0;
img = zeros( geom.reconSize, 'single' );


% apply ramp filter and hilber transfrom to the entire sinogram
sino = cosineWeighting( sino, geom );
sino2 = sino;
for iview = 1:geomSlice.noViews
    
    if tiltedFilter
        proj = sino(:, :, iview );
        % tilte the sinogram
        proj = tilteProjction(  proj, vv , vt );
        %filtering
        proj1 = filterFreq( proj, H, 2 );
        proj2 = hilberTransform( proj, 2, 2 );
        % tilte the sinogram back
        sino(:, :,  iview)  = tilteProjction(  proj1, vt , vv );
        sino2(:, :,  iview) = tilteProjction(  proj2, vt , vv ) * rampFilterGain;
    else
        proj = sino(:, :, iview );
        sino(:, :,  iview)  = filterFreq( proj, H, 2 );
        sino2(:, :,  iview) = hilberTransform( proj, 2, 2 ) * rampFilterGain;
    end
    
end

% compute redundant weighted for one slice reconstruction
redundantWeights = computeHelicalWeights2D( geom, noViewsPerSlice );


fprintf( '\tSlice: ');
for iz = 1: nz
    
    if mod(iz, 10 ) == 0
        fprintf( '%i/(%i), ', iz, nz );
    end
    
    z = ( iz - (nz+1)/2 + oz ) * dz ;
    [ ~, iview_middle ] = min(abs(couchZ-z));
    
    % determine the begin and end veiw
    iview_start = iview_middle - noViewsPerSlice/2 + 1;
    iview_stop  = iview_middle + noViewsPerSlice/2;
    
    iview_start = max( iview_start, 1 );
    iview_stop  = min( iview_stop, geom.noViews );
    
    geomSlice.reconOffset(3)   = - z / dz;
    geomSlice.noViews   = iview_stop - iview_start + 1;
    geomSlice.betas     = geom.betas(   iview_start : iview_stop );
    geomSlice.couchZ    = geom.couchZ(  iview_start : iview_stop );
    geomSlice.shortScan = 1;
    
    % not enough projection
    if iview_stop - iview_start + 1 < noViewsPerSlice
        weightViewOffset = floor( ( noViewsPerSlice - geomSlice.noViews ) / 2 );
    else
        weightViewOffset = 0;
    end
    
    
    sinoSlice = zeros( [geom.detSize(2) geom.detSize(1) geomSlice.noViews], 'single' );
    
    % apply weighting to the segement of sinogram for single slice reconstruction
    for iview = 1:geomSlice.noViews
        
        weightsProj = redundantWeights(:,:,iview + weightViewOffset );
        dRedundantWeights = imfilter( weightsProj, [-1 0 1], 'same',  'replicate' ) /  2 ;
        
        sinoSlice(:, :,  iview) =  weightsProj.*  sino(:, :, iview + iview_start - 1 ) ...
            + 2 *  dRedundantWeights .* sino2(:, :, iview + iview_start - 1 );
        
    end
    
    % finally, single slice FBP recon
    img(:,:,nz-iz+1) = backProjectMex(sinoSlice, geomSlice,  1, 0, 'back,hpd' ) * dbeta * 10;
    
end

end
