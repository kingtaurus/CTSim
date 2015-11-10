function [ img, sinoNormalized, imgMetalMap ] = reconFBPNormalizedMetalArtifactReduction(...
    sino, geom, spectrum, sinoMapUnknown, addMetalBack )
%function [ img, sinoCorrected, sinoMapUnknown, imgMetalMap ] =
%   reconFBPNormalizedMetalArtifactReduction( sino, geom, spectrum,
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
% Based on:
% Meyer, Esther, Rainer Raupach, Michael Lell, Bernhard Schmidt, and Marc
% Kachelrieß. “Normalized Metal Artifact Reduction (NMAR) in Computed
% Tomography.” Medical Physics 37, no. 10 (2010): 5482. doi:10.1118/1.3484090.
%
% 
% 2013.3 Meng Wu at Stanford University

if nargin < 5
    addMetalBack = false;
end

% paramters for NMAR method
nitn    = 2;

%initial adaptive segmentation thresholds for 120 kVp
clusterThresholds = [0.1 0.35 0.5 0.7];
metalThreshold = 1;

% segment the metal out of the image
imgFBP = reconFBP( sino, geom);
imgMetalMap = imgFBP > metalThreshold;
imgMetalMap = imdilate( imgMetalMap, ones(3));

% reconstruct estimated image use FBP
img = reconInterpFBP( sino, geom, spectrum, sinoMapUnknown );

if sum( sinoMapUnknown(:) ) == 0
    fprintf('Warning: no metal detected. \n');
    return;
end

for itn = 1:nitn
    
    % compute synthetic image
    [clusterThresholds, clusterMeans] = segmentationAdaptiveThreshold( img, clusterThresholds);
    
    map        = img < clusterThresholds(1);
    img( map ) = 0;
    map        = clusterThresholds(1) < img & img < clusterThresholds(2);
    img( map ) = clusterMeans(2);
    
    % foward projections
    img = extendVoi( img, 2 );
    sinoPrior =  forwardProjectMex( img, geom );
    sinoPrior( sinoPrior < 0.5 ) = 0.5;   
    
    sinoNormalized = sino ./ sinoPrior;
    
    % interpolation
    sinoNormalized = sinogramInterpolation( sinoNormalized, sinoMapUnknown, 'linear' );

    % denormalize
    sinoNormalized = sinoNormalized .* sinoPrior;
      
    sinoNormalized = beamHardeningWarterCorrection(sinoNormalized, spectrum );
    
    %FBP reconstruction
    img = reconFBP( sinoNormalized, geom);
       
end

if addMetalBack
    img( imgMetalMap ) = max( img( imgMetalMap ) , imgFBP( imgMetalMap ));
end

end

