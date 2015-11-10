function I = imdilate3( I, h )


if ismatrix(I)
    I = imdilate( I, h );
elseif ndims( I ) == 3
    for iz = 1: size( I, 3)
       I(:,:,iz) = imdilate(I(:,:,iz), h);
    end

end