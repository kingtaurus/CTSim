function [img, FOV] = removeOutsideFOV( img, geom, FOV )
% function img = removeOutsideFOV( img, geom )
% 
% remove image outside the FOV
%
% Meng Wu
% 2013


x = ( (-(geom.reconSize(1)-1)/2 : (geom.reconSize(1)-1)/2 ) + geom.reconOffset(1) ) * geom.reconSpacing(1);
y = ( (-(geom.reconSize(2)-1)/2 : (geom.reconSize(2)-1)/2 ) + geom.reconOffset(2) ) * geom.reconSpacing(2);
[xx, yy] = meshgrid( x, y);


if nargin < 3
    FOV = 2 * max( abs( x ) ) * geom.SAD / geom.SDD ;
end


fovMap = ( ( xx.^2 + yy.^2 ) < (FOV/2)^2 );

for iz = 1 : size(img, 3)
    
    slice = img(:,:,iz);
    slice(~fovMap) = 0;
    img(:,:,iz) = slice;
    
end


end