function img = reconATRACT1( sino, geom, window, crop, collimationBoundaries )
% function img = reconATRACT2( sino, geom, window, crop )
%
% Input:
%
% Output:
%
% Approximated Trucation Robust Algorithm for Computed Tomography 2D
% version
%   1. Work for flat detector axial cone beam CT
%
% Based on the Xia et al. "Towards clinial application of a Laplace
% operator-based region of interest reconstruction algorithm in C-Arm CT,"
% IEEE TMI, 2014
%
% Implemented by Meng Wu at 2014.3

if nargin < 4
    if nargin < 3
        window = 'hamming';
    end
    crop = 1;
end

if nargin < 5
    collimationBoundaries = [1 geom.detSize(1) 1 geom.detSize(2)];
end

edgeWidth = 5;

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
    H       = designResigualFilter2D( geom.detSize, geom.detSpacing, window, crop );
    gamma   = atan( - u / SDD );
    Tau     = abs(geom.betas(end) - geom.betas(1)) - pi;
else
    error('Only works for flat panel detecotor for now!\n');
end

% Hramp = designFilter2(du, window, nu, 1);
% Hlap = fft( [-1 2 -1], length(H) );
% Heff = Hlap .* H;
% 
% close all;
% figure; plot( Hramp ); title 'ramp';
% figure; plot( abs( Hlap ) );  title 'laplace';
% figure; plot( real( fft( H ) ) ); title 'residue';
% figure; plot(  abs( Heff ) ); title 'eff';


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
    
    proj = LaplaceFiltering1D( proj, collimationBoundaries, edgeWidth );
    proj = filterFreq( proj, H, 2 ) ;
    
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

function viewFlit = LaplaceFiltering1D( view, collimationBoundaries, edgeWidth  )

h = [-1 2 -1];
%h = [0.5 1 0.5; 1 -6 1; 0.5 1 0.5];

viewFlit = imfilter( view, h,  'replicate', 'same' );

viewFlit(:,1:collimationBoundaries(1)+edgeWidth) = 0;
viewFlit(:,collimationBoundaries(2)-edgeWidth+1:end) = 0;
viewFlit(1:collimationBoundaries(3)+edgeWidth, :) = 0;
viewFlit(collimationBoundaries(4)-edgeWidth+1:end, :) = 0;

end

function H = designResigualFilter2D( detSize, detSpacing, window, crop )

pnu = max(64,2^nextpow2(2*detSize(1)-1));

wu =  ( -(pnu-1)/2:(pnu-1)/2 ) / (2 * detSpacing(1)  * detSize(1)  );

H = fftshift( abs( wu ) ./ ( wu.^2 ) );
W = designWhindow2( window, detSize(1), crop) ;

H = H .* W';

% for iv = 1 : pnv
%     H(iv,:) =   H(iv,:) .* W';
% end


end