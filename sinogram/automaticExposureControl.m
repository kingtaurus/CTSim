function [ effectivePhotons, tubeCurrentProfile ] = automaticExposureControl( effectivePhotons, spectrum, desiredPhotonCounts, unattenPhotonCounts )
% function [ effectivePhotons, tubeCurrentProfile ] = automaticExposureControl( effectivePhotons, spectrum, desiredPhotonCounts, unattenPhotonCounts )
%
% Based on Gies et al, Medical Physics, 26, 2235 (1999)
%
% Meng Wu at Stanford University
% 2014


referenceType = 'min_column';
%correctType = 'next';
correctType = 'current';

% maximum and minimum tube output constrains
maximumScaleFactor = spectrum.maxPhotonsPerPixel / spectrum.photonsTotalOriginal;
minimumScaleFactor = 0.1;

if ndims(effectivePhotons) == 3
    nu = size( effectivePhotons, 2);
    noViews = size( effectivePhotons, 3 );
else
    nu = size( effectivePhotons, 1);
    noViews = size( effectivePhotons, 2 );
end

iu = round( nu / 2 );

if spectrum.automaticExposureControl == 0
    
    tubeCurrentProfile = ones( noViews, 1 );
    fprintf('\tDo not use automatic exposure control.\n');
    
else
    
    fprintf('\tUse automatic exposure control using %f power of attenuation. \n', spectrum.automaticExposureControl);
    fprintf('\t\tExpected %g photon counts at 25 cm soft tissue thickness.\n', desiredPhotonCounts);
    fprintf('\t\tUse %s as reference area.\n', referenceType);
    
    tubeCurrentProfile = ones( noViews, 1 );
    
    % calculate base attenuation
    baseAtten = log( unattenPhotonCounts / desiredPhotonCounts );
    
    for iv = 1 : noViews
        
        if ndims(effectivePhotons) == 3
            
            if strcmpi( referenceType, 'central_column' ) % use central ray
                
                referenceArea = effectivePhotons(:, iu-3:iu+4, iv );
                referenceAtten = log(  unattenPhotonCounts / mean( referenceArea(:) ) );
                
            elseif strcmpi( referenceType, 'min_column' ) % use column with minimum counts
                
                referenceArea = min( mean( effectivePhotons(:, :, iv ) ) );
                referenceAtten = log(  unattenPhotonCounts / mean( referenceArea ) );
                
            elseif strcmpi( referenceType, 'central_tile' ) % use central 50 X 50 pixels
                
                referenceArea = mean( effectivePhotons( end/2-24:end/2+25, end/2-24:end/2+25, iv ) );
                referenceAtten = log(  unattenPhotonCounts / mean( referenceArea(:) ) );
                
            else
                
                fprintf('Warning: unknow AEC referenceType. \n');
                referenceAtten = baseAtten;
                
            end
            
            scaleFactor = ( referenceAtten / baseAtten )^spectrum.automaticExposureControl;
            % maximum and minimum tube output constrains
            scaleFactor = min(  maximumScaleFactor, scaleFactor );
            scaleFactor = max(  minimumScaleFactor, scaleFactor );
            
            if strcmpi( correctType, 'next') && iv < noViews
                
                % scale the effective Photons in next frame
                tubeCurrentProfile( iv + 1 ) = scaleFactor;
                effectivePhotons(:,:,iv + 1) = scaleFactor * effectivePhotons(:,:, iv + 1) ;
                
            else
                % scale the effective Photons in current frame (ideal case)
                tubeCurrentProfile( iv ) = scaleFactor;
                effectivePhotons(:,:,iv ) = scaleFactor * effectivePhotons(:,:, iv ) ;
                
            end
            
            
        else
            
            referenceArea = effectivePhotons( iu-8:iu+8, iv );
            referenceAtten = log(  unattenPhotonCounts / mean( referenceArea(:) ) );
            
            scaleFactor = ( referenceAtten / baseAtten )^spectrum.automaticExposureControl;
            
            % maximum and minimum tube output constrains
            scaleFactor = min(  maximumScaleFactor, scaleFactor );
            scaleFactor = max(  minimumScaleFactor, scaleFactor );
            
            if strcmpi( correctType, 'next') && iv < noViews
                % scale the effective Photons in next frame
                tubeCurrentProfile( iv + 1 ) = scaleFactor;
                effectivePhotons(:, iv + 1 ) = scaleFactor * effectivePhotons(:, iv + 1 );
            else
                tubeCurrentProfile( iv ) = scaleFactor;
                effectivePhotons(:, iv ) = scaleFactor * effectivePhotons(:, iv ) ;
            end
            
        end
        
    end
    
end

