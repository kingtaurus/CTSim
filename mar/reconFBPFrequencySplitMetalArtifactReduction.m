function [ img, sinoCorrected, imgMetalMap ] = reconFBPFrequencySplitMetalArtifactReduction(...
    sino, geom, spectrum, sinoMapUnknown, sigma  )
%function [ img, sinoCorrected, sinoMapUnknown, imgMetalMap ] =
%   reconFBPFrequencySplitMetalArtifactReduction( sino, geom, spectrum,
%   sinoMapUnknown, addMetalBack )  
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
% Based on: Meyer, Esther, Rainer Raupach, Michael Lell, Bernhard Schmidt,
% and Marc Kachelrieß. “Frequency Split Metal Artifact Reduction (FSMAR) in
% Computed Tomography.” Medical Physics 39, no. 4 (April 2012): 1904–16.
% doi:10.1118/1.3691902.
% 
% 2013.3 Meng Wu at Stanford University

if nargin < 5
    sigma = 30 / 2.355; % 3 /cm FWHM of Guassian
end

[ imgNMAR, sinoCorrected, imgMetalMap ] = reconFBPSegmentationMetalArtifactReduction(...
    sino, geom, spectrum, sinoMapUnknown, 1 );

h = fspecial('gaussian', 16, 2 * pi / ( sigma * geom.reconSpacing(1) ) );

% frequency splits of the NMAR image
imgNMARLow  = imgNMAR;
for iz = 1:size( imgNMAR, 3 )
    imgNMARLow(:,:,iz) =  imfilter( imgNMAR(:,:,iz), h, 'replicate');
end
imgNMARHigh = imgNMAR - imgNMARLow;

% frequency splits of the FBP image
imgFBPHigh  = reconFBP( sino, geom);
for iz = 1:size( imgNMAR, 3 )
    imgFBPHigh(:,:,iz) = imgFBPHigh(:,:,iz) -  imfilter( imgFBPHigh(:,:,iz), h, 'replicate');
end

weights = single( imgMetalMap ) ;
h2 = fspecial('gaussian', 25,   10 * pi / ( sigma * geom.reconSpacing(1) ) );
for iz = 1:size( imgNMAR, 3 )
    weights(:,:,iz) =  imfilter( weights(:,:,iz), h2, 'replicate');
end

% combine frequency splitted images
img = weights .* imgFBPHigh +  ( 1 - weights ) .* imgNMARHigh + imgNMARLow;
