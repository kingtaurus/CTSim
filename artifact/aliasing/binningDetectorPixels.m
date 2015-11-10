function [sinoOut, geomOut] = binningDetectorPixels( sinoIn, geomIn, binSize )
% Combine multiples rows and columns of the detector for FDK method
% Naively conv with a rectangular function and then downsample it.
% Careful there will be aliasing !
%
% Meng Wu at Stanford Univeristy
% 2012 - 2013

if mod( geomIn.detSize(1), binSize(1) ) ~= 0
    fprintf( 'Warning: the detector size has to be a multiple of binning ratio! \n');
end


geomOut = geomIn;

if length( size( sinoIn )) == 2
    
    geomOut.detSpacing(1)   = geomIn.detSpacing(1) * binSize(1);
    geomOut.detSize(1)      = geomIn.detSize(1) / binSize(1);
    geomOut.detOffset(1)    = geomIn.detOffset(1) / binSize(1);
    
    sinoOut = detectorPixelBinning( sinoIn, [binSize(1) 1] );

elseif length( size( sinoIn )) == 3
    
    geomOut.detSpacing(1)   = geomIn.detSpacing(1) * binSize(1);
    geomOut.detSpacing(2)   = geomIn.detSpacing(2) * binSize(2);
    geomOut.detSize(1) = geomIn.detSize(1) / binSize(1);
    geomOut.detSize(2) = geomIn.detSize(2) / binSize(2);
    geomOut.detOffset(1)    = geomIn.detOffset(1) / binSize(1);
    geomOut.detOffset(2)    = geomIn.detOffset(2) / binSize(2);
    
    sinoOut = zeros( [ geomOut.detSize(2)  geomOut.detSize(1) size(sinoIn, 3) ], 'single');
    for iv = 1: size( sinoIn, 3)
        sinoOut(:,:,iv) = detectorPixelBinning( sinoIn(:,:,iv), binSize );
    end
    
end

end
