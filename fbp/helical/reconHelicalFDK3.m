function img = reconHelicalFDK3( sino, geom, window, crop, mode )
% function [img, sino] = reconHelicalFDK( sino, geom, window, crop, mode )
%   
% Not tilted filtering in this version. Very naive fdk  
% FDK helical CT recontruction not efficient in weighting and filtering, I
% guess it more robust in noisy data than the other version
%
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
%
% Meng Wu, Stanford University, 2014

if nargin < 5
    if nargin < 4
        if nargin < 3
            window = 'hamming';
        end
        crop = 1;
    end
    mode = 0;
end

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
SAD         = geom.SAD;
noViewsTurn = geom.noViewsTurn;
noViews     = geom.noViews;
h           = geom.couchSpacing;
couchZ      = geom.couchZ;

% array of detector positions in mm
u = ( ( -(nu-1)/2:(nu-1)/2) + ou ) * du;
v = ( ( -(nv-1)/2:(nv-1)/2) + ov ) * dv ;
[uu, vv ] = meshgrid( u, v);


% weighting factor and ramp filter
if geom.flatPanel
    weight  = SAD * sqrt( 1.0 + ( vv.^2 / SDD^2 ) ) ./ sqrt(  uu.^2 + vv.^2 + SDD^2 );
    H       = du*designFilter2(du, window, nu, crop);
    gamma   = atan( - u / SDD );
else
    weight  = SAD / SDD * cos( uu ./ ( SDD * sqrt( 1 + (vv / SDD ).^2 )) );
    H       = du*designEquiangularFilter2(du, SDD, window, nu, crop);
    gamma   = - u / SDD ;
end
WH = designWhindow2( window, nu, crop);
rampFilterGain = sum( real( fft(H)));

dbeta = abs(geom.betas(end) - geom.betas(1)) / (noViews-1);
vt   = vv - repmat( gamma , nv, 1) / dbeta * h;

if mode == 0 % pi + fan segments
    noViewsPerSlice = round( noViewsTurn /4  +  max(abs(gamma))/dbeta + 2 ) * 2;
elseif mode == 1 % 1.5 pi segments
    noViewsPerSlice = round( 3 * noViewsTurn / 4 ) ;
elseif mode == 2 % 3 pi segments
    noViewsPerSlice = round( 3 * noViewsTurn / 2 ) ;
end

geomSlice = geom;
geomSlice.reconSize(3)  = 1;
geomSlice.helicalScan   = 0;
img = zeros( geom.reconSize, 'single' );

for iz = 1: nz
    
    % slice position to reconstruct
    z = ( iz - (nz+1)/2 + oz ) * dz ;
    [ ~, iview_middle ] = min(abs(couchZ-z));

    % determine the begin and end veiw
    iview_start = find( couchZ > z + min(v) / dbeta * h, 1, 'first' );
    iview_stop  = find( couchZ < z + max(v) / dbeta * h, 1, 'last'  );
    
    iview_start = max( iview_start, 1 );
    iview_stop  = min( iview_stop, geom.noViews );
    
    iview_start = max( iview_start, iview_middle - noViewsPerSlice/2 + 1 );
    iview_stop  = min( iview_stop,  iview_middle + noViewsPerSlice/2 );
    
    % not enough projection to do a short scan so just scape this slice
    if iview_stop - iview_start < geom.noViewsTurn /2
        continue;
    end
    
    geomSlice.reconOffset(3)   = - z / dz ;
    geomSlice.noViews   = iview_stop - iview_start + 1;
    geomSlice.betas     = geom.betas(   iview_start : iview_stop );
    geomSlice.couchZ    = geom.couchZ(  iview_start : iview_stop );
    geomSlice.shortScan = 1;
    
    % parameters for short scan weights
    Tau                 = abs(geomSlice.betas(end) - geomSlice.betas(1)) - pi;
    betaMin             = min(geomSlice.betas);
    
    sinoSlice = zeros( [nv nu geomSlice.noViews], 'single' );
    %linear interpolate sinogram for a slice
    for iview = 1:geomSlice.noViews
        
        weightShortScan =  shortScanSliverWeight( gamma, geomSlice.betas(iview) - betaMin,  Tau );
        proj = repmat( weightShortScan, nv, 1) .* weight .* sino(:, :, iview + iview_start - 1 );
        
        for iv = 1 : nv
            projs = squeeze( proj(iv,:) ) ;
            % filter projection slice
            proj( iv, :) = filterFreq( projs, H ) ;
        end
        sinoSlice(:, :,  iview) = proj;
        
    end
    
    % finally, single slice FBP recon
    img(:,:,nz-iz+1) = backProjectMex(sinoSlice, geomSlice,  1, 0, 'back,pd' ) * dbeta * 10;
end

end


