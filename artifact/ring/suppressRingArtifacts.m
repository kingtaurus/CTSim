function img = suppressRingArtifacts( img, geom, threshold, radius )
% Suppress ring artifacts in the CT volume
%
% Meng Wu at Stanford University
% 2013


if nargin < 4
    radius = geom.reconSize(1) * geom.reconSpacing(1) / 2;
end

gf = fspecial('gaussian',[5 5], 0.6);
for iz = 1 : size( img, 3)
    imgRing = getRingArtifactPolar( squeeze(img(:,:,iz)), geom, threshold, radius );
    
    % low pass filter the ring artifacts
    imgRing = imfilter(imgRing, gf);
    img(:,:,iz) = img(:,:,iz) - imgRing;
end




end