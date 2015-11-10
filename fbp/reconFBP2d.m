function img = reconFBP2d( sino, geom, window, crop )
% Fan-beam CT filtered backpeoject reconstruction
% inout:
%       sino    - log sinogram
%       geom    - geometry parameters
%       window  - window function
% output:
%       img     - reconstructed attenuation image (1/cm)
%
% Lei Zhu, Ph.D., Stanford University (copyright reserved, May 1, 2007)
% modifications by Andreas Keil, Stanford University, 2010-11-21 - 2011-11-28
% modifications by Meng Wu, Stanford University, 2012-9

if nargin < 4
    if nargin < 3
        window = 'hamming';
    end
    crop = 1;
end


nu = geom.detSize;
du = geom.detSpacing;
noViews = geom.noViews;
SDD = geom.SDD;
SAD = geom.SAD;
minbetas = min( geom.betas );

% array of detector positions in mm
u = (( -(nu-1)/2:(nu-1)/2)+ geom.detOffset ) * du  ;

% cosine weighting factors and filter
if geom.flatPanel
    weight  = SAD ./ sqrt( u.^2 + SDD^2 );
    H       = du*designFilter2(du, window, nu, crop);
    gamma   = atan( - u / SDD );
    Tau     = abs( geom.betas(end) - geom.betas(1) ) - pi;
else
    weight  = SAD / SDD * cos( u / SDD );
    H       = du*designEquiangularFilter2(du, SDD, window, nu, crop);
    gamma   = - u / SDD ;
    Tau     = abs( geom.betas(end) - geom.betas(1) ) - pi;
end

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
    sino(:, view) = filterFreq( proj, H ) ;
    
end

dbeta = abs(geom.betas(end) - geom.betas(1)) / (noViews-1);
if ~geom.shortScan
    dbeta = dbeta / 2;
end

%Scale reconstructed attenuation image to (1/cm)
img = backProjectMex(sino, geom,  1, 0, 'back,pd' ) * dbeta * 10;


img( ~ geom.map ) = 0;
end