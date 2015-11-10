function img = reconShiftInvariantFBPShortScan( sino, geom, window, crop )
% Cone-beam CT filtered backpeoject reconstruction using FDK method
% input:
%       sino    - log sinogram
%       geom    - geometry parameters
%       window  - window function
% output:
%       img     - reconstructed attenuation image (1/cm)
%
% Based on Zhu et al "A short-scan reconstruction for cone-beam CT using
% shift-invariant FBP and equal weighting" MedPhys, 2007
%
% Meng Wu, Stanford University, 2014-3

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

% array of detector positions in mm
u = ( ( -(nu-1)/2:(nu-1)/2) + ou ) * du ;
v = ( ( -(nv-1)/2:(nv-1)/2) + ov ) * dv ;
[uu, vv ] = meshgrid( u, v);

% weighting factor 1
weight = SAD * sqrt( 1 + ( vv.^2 / SDD^2 ) ) ./ sqrt( uu.^2 + vv.^2 + SDD^2 );

% filter and window design
if geom.flatPanel
    H = du * designFilter2(du, window, nu, crop);
else
    H = du*designEquiangularFilter2( du, SDD, window, nu, crop);
end
WH = designWhindow2( window, nu, crop);

rampFilterGain = sum( real( fft(H)));
betaMiddle  = mean( geom.betas );

sino = binningDetecotRows( sino, round(  geom.reconSpacing(3) / geom.detSpacing(2)  ) );
sino2 = sino;

% filter projections
for view = 1:noViews
    % get current view's sinogram
    proj = sino(:, :, view);
    
    % cosine weighting
    proj = weight .* proj;
    proj2 = proj;

    % ramp filtering and differential
    for iv = 1 : nv
        projs = squeeze( proj(iv,:) ) ;
        proj( iv, :) = filterFreq( projs, H ) ;
        proj2( iv, :) = filterFreq( projs, WH ) ;
    end
    
    % get current view's sinogram
    sino(:, :,  view)   = proj;
    sino2(:, :,  view) = imfilter( proj2, [-1 0 1], 'same',  'replicate') /  2 ;
end

%Scale reconstructed attenuation image to (1/cm)
img1 = backProjectMex(sino, geom,  1, 0, 'back,pd' );

geom2 = geom;
%geom2.reconSize(1:2) = geom.reconSize(1:2) ;
geom2.betas = geom2.betas - betaMiddle;

%back project a rotated differential image
img2 = backProjectMex(sino2, geom2,  1, 0, 'back,pd' );

% hilbert transform
for iz = 1 : size( img2, 3 )
    
    temp = img2(:,:,iz);
    for j = 1 :  size( temp, 2)
        % zero padding before hilbert transform
        temp(:,j)= hilberTransform( temp(:,j), 2 );
    end

    % rotate back
    img2(:,:,iz) = imrotate( temp, - betaMiddle * 180 / pi, 'crop', 'cubic');
end

% combine images from two processes
img = 0.5 * ( img1 + img2 * rampFilterGain)  * abs(geom.betas(end) - geom.betas(1)) / noViews * 10 ;


