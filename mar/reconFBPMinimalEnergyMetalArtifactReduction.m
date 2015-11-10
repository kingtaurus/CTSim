function [ img, sinoCorrected, imgMetalMap ] = reconFBPMinimalEnergyMetalArtifactReduction(...
    sino, geom, spectrum, sinoMapUnknown, addMetalBack )
% function [ img, sinoCorrected, imgMetalMap ] = reconFBPMinimalEnergyMetalArtifactReduction(...
%    sino, geom, spectrum, sinoMapUnknown, addMetalBack )
% FBP recosntruction using normalized metal artifact reduction algorithms
% with second pass correction
% input:
%       sino            - attenuation sinogram
%       geom            - system geometry
%       spectrum        - X-ray spectrum infomation of KeV
%       sinoMapUnknown  - map of metal-affected pixels
%       addMetalBack    - add metal image back to reconstruction (default false)
% output:
%       img             - reconstruction result
%       sinoCorrected   - NMAR corrected sinogram 
%       imgMetalMap     - map of metal object in image
% 
% Based on: 
% Koehler, Thomas, Bernhard Brendel, and Kevin M. Brown. “A New Method for
% Metal Artifact Reduction in CT.” In The Second International Conference
% on Image Formation in X-Ray Computed Tomography, 29–32, 2012. 
% 
%   This is the best semi-iterative method so far... 
%
% 2013.3 Meng Wu at Stanford University

if nargin < 5
   addMetalBack = false; 
end

[ img, sinoCorrected, imgMetalMap ] = reconFBPSegmentationMetalArtifactReduction(...
    sino, geom, spectrum, sinoMapUnknown, false );

img = extendVoi( img, 2 );
img( imgMetalMap ) =  img( imgMetalMap );

% post metal artifacts reduction by minimizing local roughness
sinoDiff = sino - forwardProjectMex( img, geom );
gf = fspecial('gaussian',[3 3], 0.5);
for iv = 1 : size( sinoDiff, 3 )
    sinoView = imfilter(sinoDiff(:,:,iv), gf, 'replicate');
    sinoDiff(:,:,iv) = sinoView .* imdilate( sinoMapUnknown(:,:,iv), ones(3) );
end

correctionMap = 0.05 < img & img < 0.35;
imgError = reconFBP( sinoDiff , geom);

for iz = 1 : size( imgError, 3)
    map = imerode( correctionMap(:,:,iz), ones(5) );
    img(:,:,iz) = minimizeLocalStd( img(:,:,iz), imgError(:,:,iz), 11, map );
    %img(:,:,iz) = minimizeLocalEntropy( img(:,:,iz), imgError(:,:,iz), 11, map );
end


if addMetalBack
    metalThreshold = 3;
    % segment the metal out of the image
    imgFBP = reconFBP( sino, geom);
    imgMetalMap = imgFBP > metalThreshold;
    imgMetalMap = imdilate( imgMetalMap, ones(3));
    img( imgMetalMap ) = max( img( imgMetalMap ) , imgFBP( imgMetalMap ));
end


end


