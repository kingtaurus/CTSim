function imgOut = imageFilter2D( imgIn, h )

imgOut = zeros( size(imgIn), class(imgIn) );

for iz = 1 : size( imgIn, 3 )
    
    imgOut(:,:,iz) = imfilter( squeeze( imgIn(:,:,iz) ), h, 'replicate',  'same' );
    
end



end