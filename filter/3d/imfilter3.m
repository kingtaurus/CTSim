function I = imfilter3( I, h)

if ismatrix(I)
    I = imfilter( I, h );
elseif ndims( I ) == 3
    for iz = 1: size( I, 3)
       I(:,:,iz) = imfilter(I(:,:,iz), h, 'same' );
    end

end