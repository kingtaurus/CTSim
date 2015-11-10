function countingNoise = addDetectorNoisePowerSpectrum( countingNoise, geom  )
%% add countingNoise power spectrum


if isfield( geom, 'detNPS' ) && geom.detNPS > 0
    % FWHM to gaussian blurr
    
    if ismatrix( countingNoise )
        sigma = geom.detNPS / 2.355 / geom.detSpacing(1) ;
        width = ceil( sigma * 3 ) * 2 + 1;
        h = fspecial('gaussian', [width 1], sigma );
        
        countingNoise = imfilter( countingNoise, h, 'replicate', 'same' );
        
    else
        
        sigma1 = geom.detNPS  / 2.355 / geom.detSpacing(1) ;
        sigma2 = geom.detNPS  / 2.355 / geom.detSpacing(2) ;
        width1 = ceil( sigma1 * 3 ) * 2 + 1;
        width2 = ceil( sigma2 * 3 ) * 2 + 1;
        
        h1 = fspecial('gaussian', [1 width1], sigma1 );
        h2 = fspecial('gaussian', [width2 1], sigma2 );
        
        h = h2 * h1;
        
        for iv = 1:size(countingNoise, 3)
            countingNoise(:,:,iv) = imfilter( countingNoise(:,:,iv), h, 'replicate', 'same' );
        end
        
    end
    
    
    
    
end