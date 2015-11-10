function sinoMetalMap = findUntrustedPixel( sino, geom, spectrum, threshold, imgFBP )
% Find the metal pixels in the singoram and return binary maps
%   input:
%       sino        - sinogram
%       geom        - geometry parameters
%       spectrum    - spectrum parameters for compute threshold
%       threshold (default 1000)
%       sinogramBased - 0/1
%   output:
%       sinoMetalMap - sinogram metal map
%       imgMetalMap  - image metal map
%
% Meng Wu @stanford
% 2012-2013


if nargin < 4
    threshold = 500;
end

metalMarginPixels = 3;
metalMarginAngles = 2;

fprintf( 'Find the metal pixels ' );

% chech if it is attenuation coefficient line integral sinogram
if max( sino(:) ) > 1000
    sino = computeSinogramAttunation( sino, spectrum );
end

if threshold > 100
    
    
    fprintf( 'in singoram space... \n' );
    
    % determine threshold for untrusted keV attenuation integral values
    spetrumThreshold = threshold * ...
        spectrum.photonsPerEnergyBin / sum( spectrum.photonsPerEnergyBin);
    
    intensityThreshold = sum(hDetectorEnergyResponseFunction( ...
        spectrum.DQE * spetrumThreshold, spectrum.energyBinLabels));
    
    sourceIntensity = sum( hDetectorEnergyResponseFunction( ...
        spectrum.DQE * spectrum.photonsPerEnergyBin, spectrum.energyBinLabels));
    
    intensityThresholdRatio = intensityThreshold / sourceIntensity;
    
    AttThreshold = -log(intensityThresholdRatio);
    
    
    fprintf('\tNot trusting transmitted keV photon numbers less than %g, \n\t i.e. keV attenuation values greater than %g.\n', ...
        threshold, AttThreshold);
    
    
    sinoMetalMap = sino > AttThreshold;
    
else
    
    metalThreshold = threshold;
    
    fprintf( 'in image space with threshold %.2f. \n' , metalThreshold);
    
    
    if nargin < 5
        sino = medianFilterSino( sino, 3 );
        imgFBP = reconFBP( sino, geom, 'ram-lak');
    end
    
    %imgFBP = extendVoi( imgFBP, 1 );
    imgMetalMap = imgFBP > metalThreshold;
    
    sinoMetal = forwardProjectMex( single(imgMetalMap), geom );
    sinoMetalMap = logical( sinoMetal > 0.1 );
    
end

%% dilate metal map

sinoMetalMapTemp = sinoMetalMap;

if length( geom.reconSize ) == 3 && length( geom.detSize ) == 2
    
    se =  ones( metalMarginPixels );
    for iv = 1 : geom.noViews
        
        ivp = iv - 1;
        if ivp == 0 && ~ geom.shortScan
            ivp = geom.noViews;
        else
            ivp = iv;
        end
        
        ivn = iv + 1;
        if ivn > geom.noViews && ~ geom.shortScan
            ivn = ivn - geom.noViews;
        else
            ivn = iv;
        end
        
        if metalMarginAngles > 1
            
            sinoMetalMap(:,:,iv) = logical(imdilate(sinoMetalMapTemp(:,:,iv), se)) | ...
                logical(imdilate(sinoMetalMapTemp(:,:,ivp), se)) | logical(imdilate(sinoMetalMapTemp(:,:,ivn), se));
        else
            
            sinoMetalMap(:,:,iv) = logical(imdilate(sinoMetalMapTemp(:,:,iv), se));
        end
        
        
    end
    
    
elseif length( geom.reconSize ) == 2 && length( geom.detSize ) == 1
    
    se = ones( metalMarginPixels,  metalMarginAngles);
    sinoMetalMap = logical(imdilate(sinoMetalMap, se));
    
else
    error('Wrong geometry!\n');
end

noUntrustedPixels = sum(sinoMetalMap(:));
ratioUntrustedPixels = noUntrustedPixels / numel(sinoMetalMap);

fprintf('\tThis affects %g pixels or %.3f of the total keV input data.\n',  noUntrustedPixels, ratioUntrustedPixels);
fprintf('Done! \n\n')


end