% function img = reconATRACT2( sino, geom, window, crop )
% Approximated Trucation Robust Algorithm for Computed Tomography Meng's
% version
% Input:
%
% Output:
%
% Notes:
%   1. Work for flat detector axial cone beam CT
%   2. Instead of using a Laplace operator and a carefully designed the
%   residue filter, I use a derivative operator plus hilbert transform. The
%   result is quite similar to the original ATRACT algorithm. 
%   3. The real magic is when apply the local numerical filter, the edge is 
%   replaced by a small number.
%   4. I guess this would have a effect to the image uniformity. The is 
%
% Based on the Xia et al. "Towards clinial application of a Laplace
% operator-based region of interest reconstruction algorithm in C-Arm CT,"
% IEEE TMI, 2014
%
%
% Implemented by Meng Wu at 2014.3

function img = reconATRACTm( sino, geom, window, crop, collimationBoundaries )

if nargin < 4
    if nargin < 3
        window = 'hamming';
    end
    crop = 1;
end

if nargin < 5
    collimationBoundaries = [1 geom.detSize(1) 1 geom.detSize(2)];
end

edgeWidth = 3;

nu = double( geom.detSize(1) );
nv = double( geom.detSize(2) );
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
    weight  = SAD * sqrt( 1.0 + ( vv.^2 / SDD^2 ) ) ./ sqrt(  uu.^2 + vv.^2 + SDD^2 );
    W = designWhindow2( window, detSize(1), crop) ;
    gamma   = atan( - u / SDD );
    Tau     = abs(geom.betas(end) - geom.betas(1)) - pi;
else
    error('Only works for flat panel detecotor for now!\n');
end

% filter projections
for view = 1:noViews
    
    % get current view's sinogram and cosine weighting
    proj = weight .* sino(:, :, view);
    
    if geom.shortScan
        weightShortScan =  shortScanSliverWeight( gamma, geom.betas(view) - betaMin,  Tau);
        %weightShortScan = shortScanParkerWeight( gamma, geom.betas(view) - betaMin, Tau );
        %weightShortScan = shortScanModifiedParkerWeight( gamma, geom.betas(view) - betaMin, Tau );
        for iv = 1 : nv
            proj( iv, :)  = squeeze( proj(iv,:) ).* weightShortScan;
        end
    end
    
    proj = filterFreq2D( proj, W, 2 );
    [proj, edge] = firstDerivative( proj, collimationBoundaries, edgeWidth );
    proj = hilberTransform( proj + edge, 2, 2 ) ;
    % get current view's sinogram
    sino(:, :,  view) = proj;
    
end

dbeta = abs(geom.betas(end) - geom.betas(1)) / (noViews-1);

if ~geom.shortScan
    dbeta = dbeta / 2;
end

%Scale reconstructed attenuation image to (1/cm)
img = backProjectMex(sino, geom,  1, 0, 'back,pd' ) * dbeta * 10;

end

function [ viewFlit, viewEdge ] = firstDerivative( view, collimationBoundaries, edgeWidth  )

h = [-1 0 1];
%h = [0.5 1 0.5; 1 -6 1; 0.5 1 0.5];
viewFlit = imfilter( view, h,  'replicate', 'same' );

viewEdge = viewFlit;
viewFlit(:,1:collimationBoundaries(1)+edgeWidth) = 0;
viewFlit(:,collimationBoundaries(2)-edgeWidth+1:end) = 0;
viewEdge = viewEdge - viewFlit;

end
