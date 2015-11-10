function img = combineHighFreqStructures( img, imgResidue, geom, use_weighted )

if nargin < 4
    use_weighted = true;
end

nx = geom.reconSize(2);
ny = geom.reconSize(1);
dx = geom.reconSpacing(1);
dy = geom.reconSpacing(2);

% compute the local Nyquist freq
x = ( -(nx-1)/2:(nx-1)/2) * dx;
y = ( -(ny-1)/2:(ny-1)/2) * dy;

[xx, yy] = meshgrid(x , y);

if use_weighted
nr = sqrt( xx.^2 + yy.^2 ) + 20 ;
weights = sqrt( nr ./ 50 );
weights( weights(:) > 1) = 1;
weights( sqrt( xx.^2 + yy.^2 ) > geom.FOV  /2 ) = 0;
else
weights = single( sqrt( xx.^2 + yy.^2 ) <= geom.FOV /2 );
end

weights = imfilter( weights, fspecial('Gaussian', [8 8], 2 ), 'same' );

for iz = 1:size( img, 3)
        img(:,:,iz) = img(:,:,iz) + weights .* imgResidue(:,:,iz);
end
