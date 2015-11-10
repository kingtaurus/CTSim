function imgBlend = reconFBP2dMultiResolutionOld( sino, geom, window, nrsl, scale )
% Fan-beam CT filtered backpeoject reconstruction
% multi-resoltion piecewise linear frequency response FBP
% inout:
%       sino    - log sinogram
%       geom    - geometry parameters
%       window  - window function
%       nrsl    - number of resolution component
%       scale   - scale of frequency support ( default 0.5 )
% output:
%       img  - reconstructed attenuation image (1/cm)
%
% Note: 
%   1. This is the 2D version for proof the concept
%   2. The frequency are spilt equally(old version), the number of the
%   frequency component are usually higher.
%
% Meng Wu, Stanford University and University Erlangen-Nuremberg
% 2014


if nargin < 3
    window = 'hamming';
end

nu = geom.detSize;
du = geom.detSpacing;
noViews = geom.noViews;
SDD = geom.SDD;
SAD = geom.SAD;
minbetas = min( geom.betas );

nx = geom.reconSize(1);
ny = geom.reconSize(2);
dx = geom.reconSpacing(1);
dy = geom.reconSpacing(1);

x = ( -(nx-1)/2:(nx-1)/2) * dx;
y = ( -(ny-1)/2:(ny-1)/2) * dy;

dbeta = abs(geom.betas(end) - geom.betas(1)) / (noViews-1);

[xx, yy] = meshgrid(x , y);
nr = sqrt( xx.^2 + yy.^2 ) * SDD / SAD / du;
freqSupport = scale ./ ( nr * dbeta );


% array of detector positions in mm
u = (( -(nu-1)/2:(nu-1)/2)+ geom.detOffset ) * du  ;

% cosine weighting factors and filter
if geom.flatPanel
    weight  = SAD ./ sqrt( u.^2 + SDD^2 );
    H       = du*designFilter2(du, window, nu, 1);
    gamma   = atan( - u / SDD );
    Tau     = abs( geom.betas(end) - geom.betas(1) ) - pi;
else
    weight  = SAD / SDD * cos( u / SDD );
    H       = du*designEquiangularFilter2(du, SDD, window, nu, 1);
    gamma   = - u / SDD ;
    Tau     = abs( geom.betas(end) - geom.betas(1) ) - pi;
end

% multi-resolution filter
L = length(H) / nrsl / 2;
Hmr = zeros( length(H), nrsl );

for irsl = 1 : nrsl
    if irsl == 1
        Hmr( (irsl-1)*L+1 : irsl*L, irsl ) =  (L:-1:1)/L;
    else
        Hmr( (irsl-2)*L+1 : (irsl-1)*L, irsl ) =  (1:L)/L;
        Hmr( (irsl-1)*L+1 : irsl*L, irsl ) =  (L:-1:1)/L;
    end
    Hmr(:,irsl) = Hmr(:,irsl) + flipud( Hmr(:,irsl) );
    Hmr(:,irsl) = Hmr(:,irsl) .* H;
end

imgMS       = cell( nrsl, 1);
imgBlend    = zeros( geom.reconSize, 'single');

factor = 1;
if ~geom.shortScan
    factor = factor / 2;
end

sinoFilt = sino;
for irsl = 1 : nrsl
    
    H = Hmr(:,irsl);
    
    % filter projections
    for view = 1:noViews
        % get current view's sinogram
        proj = sino(:, view)';
        
        if geom.shortScan
            %short scan weighting
            weightShortScan = shortScanSliverWeight( gamma, geom.betas(view) - minbetas,  Tau);
            % cosine weighting
            proj = weight .* weightShortScan .* proj;
        else
            % cosine weighting
            proj = weight .* proj;
        end
        % get current view's sinogram
        sinoFilt(:, view) = filterFreq( proj, H ) ;
        
    end
    
    %Scale reconstructed attenuation image to (1/cm)
    img = backProjectMex(sinoFilt, geom,  1, 0, 'back,pd' ) * dbeta * 10 * factor;
    
    if irsl == 1
        freqWeighting = 1  ;
    else
        freqWeighting = 1 ./ ( 1 + exp( 6 * ( irsl / nrsl ./ freqSupport - 1  ) ) ) ;
        %freqWeighting( irsl / nrsl < freqSupport  ) = 1;
    end
    
    %imdisp( freqWeighting )
    
    imgBlend = imgBlend + freqWeighting .*  img;
    
    imgMS{irsl} = img;

end

imgBlend( ~ geom.map ) = 0;

end


