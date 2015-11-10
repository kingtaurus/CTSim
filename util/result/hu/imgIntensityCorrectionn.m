function imgAttOut =  imgIntensityCorrectionn( imgAttIn, imgGt, map  )
%   Correct image intensity based on simulated Ground truth
%   input:
%       imgAttIn - atten image to be corrected
%       imgGt    - ground truth image can be mono or polychromatoc
%       map      - map of region to have same mean intensity
%   output:
%       imgAttOut
%
% Meng Wu @ Stanford
% 2012

muTissue = mean( imgGt(map.mapTissue(:)) );


if ndims( imgAttIn ) == 3
    slice = imgAttIn(:,:,round(end/2));
    meanTissue = mean( slice( map.mapTissue(:)) );
else
    
    meanTissue = mean( imgAttIn( map.mapTissue(:)) );
end

if min( imgGt(:) ) < -500
    imgAttOut = imgAttIn * ( muTissue + 1000) / meanTissue - 1000;
else
    imgAttOut = imgAttIn * muTissue / meanTissue;
end


end