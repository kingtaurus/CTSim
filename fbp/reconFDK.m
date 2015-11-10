function img = reconFDK( sino, geom, window, crop )
% Cone-beam CT filtered backpeoject reconstruction using FDK method
% input:
%       sino    - log sinogram
%       geom    - geometry parameters
%       window  - window function
% output:
%       img     - reconstructed attenuation image (1/cm)
%
% Meng Wu, Stanford University, 2013-4

if nargin < 4
    if nargin < 3
        window = 'hamming';
    end
    crop = 1;
end

nu = geom.detSize(1);
nv = geom.detSize(2);
du = geom.detSpacing(1);
dv = geom.detSpacing(2);
ou = geom.detOffset(1);
ov = geom.detOffset(2);

noViews = geom.noViews;
SDD = geom.SDD;
SAD = geom.SAD;
betaMin = min(geom.betas);


% array of detector positions in mm
u =  ( ( -(nu-1)/2:(nu-1)/2) + ou )  * du;
v = ( ( -(nv-1)/2:(nv-1)/2) + ov ) * dv ;
[uu, vv ] = meshgrid( u, v);

% weighting factor and ramp filter
if geom.flatPanel
    weight  = SAD * sqrt( 1.0 + ( uu.^2 / SDD^2 ) ) ./ sqrt(  uu.^2 + vv.^2 + SDD^2 );
    H       = du*designFilter2( du, window, nu, crop);
    gamma   = atan( - u / SDD );
    Tau     = abs(geom.betas(end) - geom.betas(1)) - pi;  
else
    weight  = SAD / SDD * cos( uu ./ ( sqrt( SDD^2 +  vv.^2 )) );
    H       =  du*designEquiangularFilter2( du, SDD, window, nu, crop);
    gamma   = - u / SDD ;
    Tau     = abs( max( geom.betas )  - min( geom.betas )) - pi;
end

% filter projections
for view = 1:noViews
    
    % get current view's sinogram and cosine weighting
    proj = weight .* sino(:, :, view);
    
    if geom.shortScan
        %weightShortScan =  shortScanSliverWeight( gamma, geom.betas(view) - betaMin,  Tau);
        weightShortScan = shortScanParkerWeight( gamma, geom.betas(view) - betaMin, Tau );
        %weightShortScan = shortScanModifiedParkerWeight( gamma, geom.betas(view) - betaMin, Tau );
    end
    
    % ramp filtering
    for iv = 1 : nv
        
        projs = squeeze( proj(iv,:) ) ;
        if geom.shortScan
            projs = projs .* weightShortScan;
        end
        proj( iv, :) = filterFreq( projs, H ) ;
        
    end
    
    % get current view's sinogram
    sino(:, :,  view) = proj;
    
end

dbeta = abs(geom.betas(end) - geom.betas(1)) / (noViews-1);
if ~geom.shortScan
    dbeta = dbeta / 2;
end

%Scale reconstructed attenuation image to (1/cm)
img = backProjectMex(sino, geom,  1, 0, 'back,hpd' ) * dbeta * 10;

end