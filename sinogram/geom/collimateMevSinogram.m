function [sinoInterpolatedMeV, sinoMetalMap, sinoOverlapMap] = collimateMevSinogram( ...
    sinoKeV, sinoMeV, geomKeV, geomMeV, spectrumKeV, photonCountsThreshold )
% find mate
%
% Meng Wu at Stanford University
% 2012 - 2013

if nargin < 6
    photonCountsThreshold = 500;
end


collimationMarginPixels = 5;
collimationMarginAngles = 3;

% make sure keV and MeV sinograms have the same geometry
if (geomKeV.SAD ~= geomMeV.SAD) || (geomKeV.ADD ~= geomMeV.ADD)
    error('\tThis sinogram geometry conversion is not implemented yet!\n');
end

fprintf('Compute the collimated MeV sinogram ... \n');
fprintf('\tPhoton Counts Threshold = %i, Collimation Margin Pixels = %i ... \n', ...
    photonCountsThreshold, collimationMarginPixels);

sinoMetalMap = findUntrustedPixel( sinoKeV, geomKeV, spectrumKeV, photonCountsThreshold );

% sinogram map of collimated MeV acquisition
sinoCollimatedMap = sinoMetalMap;
if length( geomKeV.reconSize ) == 3 && length( geomKeV.detSize ) == 2
    
    se =  strel('disk', collimationMarginPixels);
    for iv = 1 : geomKeV.noViews
        
        ivp = iv - 1;
        if ivp == 0 && ~ geomKeV.shortScan
            ivp = geomKeV.noViews;
        else
            ivp = iv;
        end
        
        ivn = iv + 1;
        if ivn > geomKeV.noViews && ~ geomKeV.shortScan
            ivn = ivn - geomKeV.noViews;
        else
            ivn = iv;
        end
        
        sinoCollimatedMap(:,:,iv) = logical(imdilate(sinoMetalMap(:,:,iv), se)) | ...
            logical(imdilate(sinoMetalMap(:,:,ivp), se)) | logical(imdilate(sinoMetalMap(:,:,ivn), se));
    end
    
elseif length( geomKeV.reconSize ) == 2 && length( geomKeV.detSize ) == 1
    
    se = ones(collimationMarginPixels, collimationMarginAngles);
    sinoCollimatedMap = logical(imdilate(sinoMetalMap, se));
    
else
    error('Wrong geometry!\n');
end

sinoInterpolatedMeV = convertSinogramGeometry(sinoMeV, geomMeV, geomKeV );

% sinogram map of overlap area
sinoOverlapMap = logical( (~sinoMetalMap) & sinoCollimatedMap );

fprintf('done. \n');

end