function imp = cartesian2polar( imc, geom, rmin, rmax )
% Convert image from cartesian coordinate to polar coordinates
%
% Meng Wu at Stanford University
% 2013

if nargin < 3
    rmin = 0;
    rmax = max( geom.reconSize .* geom.reconSpacing ) / 2;
end

% determine sampling in azimuthal direcation
nphi = round( 2 * pi * rmax / geom.reconSpacing(1) ); 
dphi = 2 * pi / nphi;

% mesh points in polar coordinate
phi = ( 0 : dphi : (2*pi + 2*dphi) );
r   = rmin : geom.reconSpacing(1) : rmax ;
[phis, rs] = meshgrid( phi, r );

% mesh points in cartesian coordinates
xxp = cos( phis ) .* rs;
yyp = sin( phis ) .* rs;

y = ((-(geom.reconSize(1)-1)/2: (geom.reconSize(1)-1)/2) + geom.reconOffset(1)) * geom.reconSpacing(1);
x = ((-(geom.reconSize(2)-1)/2: (geom.reconSize(2)-1)/2) + geom.reconOffset(2)) * geom.reconSpacing(2);

[xxc, yyc] = meshgrid(x, y);

% 2D interpolation
imp = interp2( xxc, yyc, imc, xxp, yyp );
imp( isnan( imp) ) = 0;


end