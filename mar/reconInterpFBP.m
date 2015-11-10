function [img, sinoInterp] = reconInterpFBP( sino, geom, sinoMapUnknown, smoothMargin, addMetalBack  )
% function [img, sinoInterp, sinoMapUnknown] = reconInterpFBP( sino, geom, spectrum,
%    sinoMapUnknown, addMetalBack  )
%
% Filtered back-projection of CT recontruction using interpolation method
% to reduce metal artfacts
%
% input:
%       sino    - log sinogram
%       geom    - geometry parameters
%       spectrum - spectrum info for beam hardening correction
%       sinoMapUnknown - binary map of unknown pixels
%       addMetalBack    - add metal image back to reconstruction (default false)
% output:
%       img     - reconstructed attenuation image (1/cm)
%
% Meng Wu, Stanford University, 2013-4

if nargin < 4
    smoothMargin = 16;
end

if nargin < 5
    addMetalBack = false;
end

if sum( sinoMapUnknown(:) ) == 0
    fprintf( '\tWarning: no metal detected, use usual FBP.\n');
end

% interpolation
sinoInterp = sinogramInterpolation( sino, sinoMapUnknown, 2, smoothMargin );

img = reconFBP( sinoInterp, geom );

if addMetalBack
    metalThreshold = 1;
    % segment the metal out of the image
    imgFBP = reconFBP( sino, geom);
    imgMetalMap = imgFBP > metalThreshold;
    imgMetalMap = imdilate( imgMetalMap, ones(3));
    img( imgMetalMap ) = max( img( imgMetalMap ) , imgFBP( imgMetalMap ));
end