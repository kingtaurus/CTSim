function [sinoMapCollimatedKeV, sinoMapCollimatedMeV] = collimateSinogram(sinoMapKeV, geomKeV, geomMeV )

% sinogram map of collimated MeV acquisition
collimationMarginPixels = 1;
collimationMarginAngles = 1;
sinoMapCollimatedKeV = logical(imdilate(sinoMapKeV, ones(2*collimationMarginPixels+1, 2*collimationMarginAngles+1)));

if nargin == 3
    sinoMapCollimatedMeV = convertSinogramGeometry(sinoMapCollimatedKeV, ...
        geomKeV.noViews, geomKeV.detSize, geomKeV.detSpacing, ...
        geomMeV.noViews, geomMeV.detSize, geomMeV.detSpacing) >= 0.5;
else
    sinoMapCollimatedMeV = [];
end


end