function img = azimuthalBluring( img, geom, sigma )


rmin = 0;
rmax = max( geom.reconSize .* geom.reconSpacing ) / 2;

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


h = fspecial( 'Gaussian', [1 round(sigma * 4)], sigma );

for iz = 1 : size( img, 3)
    
    
    imp = interp2( xxc, yyc, img(:,:,iz), xxp, yyp );
    imp( isnan( imp) ) = 0;
    
    imp = imfilter( imp, h, 'circular', 'same');
    
    imc = griddata( xxp, yyp , double(imp),  xxc, yyc );
    imc( isnan( imc) ) = 0;
    
    img(:,:,iz) = imc ;
    
end

end