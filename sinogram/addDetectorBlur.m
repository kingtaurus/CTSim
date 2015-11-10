function effectivePhotons = addDetectorBlur( effectivePhotons, geom )
% function effectivePhotons = addDetectorBlur( effectivePhotons, geom )
%
% Meng Wu at Stanford University
% 2014

if isfield( geom, 'detPSF' ) && geom.detPSF > 0
    
    if length( size(effectivePhotons) ) == 2
        sigma = ( geom.detPSF + geom.focalPSF ) / 2.355 / geom.detSpacing(1) ;
        width = ceil( sigma * 3 ) * 2 + 1;
        h = fspecial('gaussian', [width 1], sigma );
        
        effectivePhotons = imfilter( effectivePhotons, h, 'replicate', 'same' );
        
    else
        
        sigma1 = ( geom.detPSF + geom.focalPSF ) / 2.355 / geom.detSpacing(1) ;
        sigma2 = ( geom.detPSF + geom.focalPSF ) / 2.355 / geom.detSpacing(2) ;
        width1 = ceil( sigma1 * 3 ) * 2 + 1;
        width2 = ceil( sigma2 * 3 ) * 2 + 1;
        
        h1 = fspecial('gaussian', [1 width1], sigma1 );
        h2 = fspecial('gaussian', [width2 1], sigma2 );
        
        h = h2 * h1;
        
        for iv = 1:size(effectivePhotons, 3)
            effectivePhotons(:,:,iv) = imfilter( effectivePhotons(:,:,iv), h, 'replicate', 'same' );
        end
        
    end
    
end

end