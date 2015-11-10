function I = imerode3( I, h )


if ismatrix(I)
    I = imerode( I, h );
elseif ndims( I ) == 3
    for iz = 1: size( I, 3)
       I(:,:,iz) = imerode(I(:,:,iz), h);
    end

end