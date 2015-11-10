function [sinoOut, geomOut] = binningDetectorColumns( sinoIn, geomIn, ratio )
% Downsample detector columns
%   Use matlab function y = decimate(x,r)
%   It lowpass filters the input to guard against aliasing and downsamples the result.
% Meng Wu at Stanford Univeristy
% 2012 - 2013


if mod( geomIn.detSize(1), ratio ) ~= 0
    error( 'The detector size has to be a multiple of binning ratio');
end


geomOut = geomIn;
geomOut.detSize(1)      = geomIn.detSize(1) / ratio; 
geomOut.detSpacing(1)   = geomIn.detSpacing(1) * ratio; 
geomOut.detOffset(1)    = geomIn.detOffset(1) / ratio;


if length( size( sinoIn )) == 2
    sinoOut = zeros( [geomOut.detSize(1) geomOut.noViews] , 'single' );
    for iview = 1: geomOut.noViews
        sinoOut(:,iview ) = decimate( double( sinoIn(:,iview) ),ratio);
    end

elseif length( size( sinoIn )) == 3
    
    sinoOut = zeros( [geomOut.detSize(2) geomOut.detSize(1) geomOut.noViews] , 'single' );
    for iview = 1 : geomOut.noViews
        for iv = 1 : geomOut.detSize(2)
            sinoOut(iv,:,iview ) = decimate( double( sinoIn(iv,:,iview) ),ratio);    
        end
    end
    
end




end
