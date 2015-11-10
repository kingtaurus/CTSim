function imgRing = getRingArtifactPolar( img, geom, threshold, radius )
% imgRing = getRingArtifactPolar( img, geom, threshold, radius )
% Get ring artifacts in the polar coordinates using median filter
%
% Meng Wu at Stanford University
% 2013


warning( 'off');
dx = geom.reconSpacing(1);

% determine the size of median filter in polar direction
Mrad = ( 11 * geom.detSpacing(1) / geom.reconSpacing(1) ) * geom.SAD / geom.SDD;
halfMrad = round(Mrad / 2);
Mrad = halfMrad* 2 + 1;

% determine the step size in polar direction where the is median filter has
% same size in azimuthal direction
Srad = round( geom.reconSpacing(1) / ( pi / 180 ) );
Srad = max( Srad, Mrad + halfMrad );

Mazi = 1;

imgRing = zeros( size(img), 'single' );

for Nrad = 1 :  Srad :  round( radius / geom.reconSpacing(1) )
    
    % convert polar coordinate
    if Nrad == 0
        imgPolar = cartesian2polar( img, geom, Nrad*dx , (Nrad + Srad + halfMrad) *dx);
    else
        imgPolar = cartesian2polar( img, geom, (Nrad - halfMrad) *dx , (Nrad + Srad + halfMrad) *dx);
    end
    
    % median filter 
    imgMedflitPolar = medfilt2( imgPolar, [ Mrad  Mazi] );
    
    imgRingArtf = imgPolar - imgMedflitPolar;
    imgRingArtf( abs(imgRingArtf) > threshold ) = 0;
    
    % convert back to cartesian coordinate
    if Nrad == 0
        imgRingArtf( end-halfMrad+1: end, :) = 0;
         imgRingCorr = polar2cartesian( imgRingArtf, geom, Nrad*dx , (Nrad + Srad + halfMrad) *dx);
    else
        imgRingArtf( end-halfMrad+1: end, :) = 0;
        imgRingArtf( 1:halfMrad, :) = 0;
         imgRingCorr = polar2cartesian( imgRingArtf, geom, (Nrad - halfMrad) *dx, (Nrad + Srad + halfMrad) *dx);
    end
    
    % get ring artifacts
    imgRing = imgRing + imgRingCorr;
    
    Mazi = Mazi + 1;
    
end

