function [ img, sinoCorrected, imgMetalMap ] = reconFBPSegmentationMetalArtifactReduction(...
    sino, geom, sinoMapUnknown, addMetalBack )
%function [ img, sinoCorrected, sinoMapUnknown, imgMetalMap ] =
%   reconFBPSegmentationMetalArtifasegCtreduction( sino, geom, spesegCtrum,
%   sinoMapUnknown, addMetalBack )
% input:
%       sino            - attenuation sinogram
%       geom            - system geometry
%       spesegCtrum        - X-ray spesegCtrum infomation of KeV
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
% No final minimization correction step.
%
% Note:
%   1. This method is very similar to the NMAR method, but gather more
%   infomation of the metal affected region via iteration. From my experience,
%   this method often out-performance the standard NMAR method.
%   2. However, this method also has its own problems: one is in case of the
%   truncation, there will be an offset between true and projected
%   sinogram. Therefore, more complicated corrections are needed.
%
% 2013.3 Meng Wu at Stanford University

if nargin < 5
    addMetalBack = false;
end

% paramters for NMAR method
alpha   = 0.8;
nitn    = 3;

%initial adaptive segmentation segThresholds for 120 kVp
segThr = [0.1 0.35 0.5 0.7];
%segThr = [0.1 0.3 0.4 ];
metalsegThreshold = 1;
%h = fspecial('Gaussian', [8 8], 4);

% segment the metal out of the image
imgFBP = reconFBP( sino, geom);
imgMetalMap = imgFBP > metalsegThreshold;
imgMetalMap = imdilate( imgMetalMap, ones(3));

% reconstruct estimated image use FBP
[img, sinoInterp] = reconInterpFBP( sino, geom, sinoMapUnknown, 0 );

if sum( sinoMapUnknown(:) ) == 0
    fprintf('Warning: no metal detected. \n');
    return;
end

for itn = 1:nitn
    
    img = immedian3( img, 3 );
    % compute synthetic image
    [segThr, segCtr] = segmentationAdaptiveThreshold( img, segThr);
    
    
    segThr(1) = min( segThr(1), segCtr(2) - 0.1 );
    segThr(2) = max( segThr(2), segCtr(2) + 0.1 );
    
    img( img < segThr(1) )      = 0;
    map = segThr(1) < img & img < segThr(2);
    map = map & ~ imgMetalMap;
    
    %map = segCtr(2) - 0.1  < img & img < segCtr(2) + 0.1;
    img( map(:) )  = softThresholding(img( map(:) ), segCtr(2), 0.05 ) ;
    %img( imgMetalMap ) = segCtr(4);
    %img = imfilter3( img, h );
    
    % foward projections
    sinoCorrected =  forwardProjectMex( img, geom );
    
    % compensation method is choosen by the user based on the truncation
    
    sinoOffset = sinogramInterpolation( sinoInterp - sinoCorrected, sinoMapUnknown, 2, 4 );
    
    % sinoOffset = sinoInterp - sinoCorrected;
    sinoCorrected( ~sinoMapUnknown )    = sinoInterp( ~sinoMapUnknown );
    sinoCorrected( sinoMapUnknown )     = sinoCorrected( sinoMapUnknown ) + sinoOffset( sinoMapUnknown );
    %sinoCorrected    = alpha * sinoCorrected + (1-alpha) * sinoInterp;
    
    %FBP reconstruction
    img = reconFBP( sinoCorrected, geom );
    
end

if addMetalBack
    img( imgMetalMap ) = max( img( imgMetalMap ) , imgFBP( imgMetalMap ));
end

end



%     for iv = 1 : size(  sinoInterp, 3 )
%
%         sinoCorrectedView   = squeeze( sinoCorrected(:,:,iv) );
%
%         sinoCorrectedView   = imfilter( sinoCorrectedView, h2 );
%         combineWeights      = imfilter( single( squeeze( sinoMapUnknown(:,:,iv) ) ), h2 );
%         sinoCorrectedView   = squeeze( sinoCorrected(:,:,iv) ) .* ( 1- combineWeights ) + sinoCorrectedView .* combineWeights;
%
%         sinoCorrected(:,:,iv) = sinoCorrectedView;
%     end


% if true
%     for iv = 1 : size(  sinoInterp, 3 )
%         sinoInterpView      = squeeze( sinoInterp(:,:,iv) );
%         sinoMapUnknownView 	= squeeze( sinoMapUnknown(:,:,iv) );
%         sinoCorrectedView   = squeeze( sinoCorrected(:,:,iv) );
%
%         %without truncation
%         offset = mean( sinoInterpView(sinoMapUnknownView(:)) - sinoCorrectedView(sinoMapUnknownView(:)) );
%         sinoCorrectedView( ~sinoMapUnknownView ) = sinoInterpView( ~sinoMapUnknownView );
%         sinoCorrectedView(sinoMapUnknownView) = alpha * ( sinoCorrectedView( sinoMapUnknownView) + offset )...
%             + (1-alpha) * sinoInterpView( sinoMapUnknownView);
%
%         sinoCorrected(:,:,iv) = sinoCorrectedView;
%     end
% else
%     %without truncation
%     offset = mean( sinoInterp(sinoMapUnknown(:)) - sinoCorrected(sinoMapUnknown(:)) );
%     sinoCorrected( ~sinoMapUnknown ) = sinoInterp( ~sinoMapUnknown );
%     sinoCorrected(sinoMapUnknown) = alpha * ( sinoCorrected( sinoMapUnknown) + offset )...
%         + (1-alpha) * sinoInterp( sinoMapUnknown);
% end


