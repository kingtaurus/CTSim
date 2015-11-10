function z = piecewiseDifferenceIso( x, imgSize,  mode )

if mode == 1
    
    Np = numel( x );
    
    if length(  imgSize ) == 2
        
        z = zeros( [ Np * 2, 1], 'single' );
        
        d = x(1:end-1, :) -  x(2:end, :);
        z( 1: numel(d) ) = d(:);
        
        d = x(:, 1:end-1) -  x(:, 2:end);
        z( Np + 1: Np + numel(d) ) = d(:);
        
    elseif length(  imgSize ) == 3
        
        z = zeros( [Np * 3, 1], 'single' );
        
        d = x(1:end-1, :, :) - x(2:end, :, :);
        z( 1: numel(d) ) = d(:);
        
        d =  x(:, 1:end-1, :) - x(:, 2:end, :);
        z( Np + 1: Np + numel(d) ) = d(:);
        
        d = x(:, :, 1:end-1) - x(:, :, 2:end);
        z( Np * 2 + 1: Np * 2 + numel(d) ) = d(:);
        
    end
    
elseif mode == 2
    
    
    z = zeros( imgSize, 'single' );
    Np = numel( z );
    
    if length(  imgSize ) == 2

        d = x( 1: numel(z(1:end-1, :)) );
        z(1:end-1, :)   = z(1:end-1, :)  + reshape( d, size(z(1:end-1, :)));
        z(2:end, :)     = z(2:end, :)  + reshape( -d, size(z(2:end, :)));
        
        d = x( Np + 1 : Np + numel(z(:, 1:end-1)) );
        z(:, 1:end-1)   = z(:, 1:end-1) + reshape( d, size(z(:, 1:end-1)));
        z(:, 2:end)     = z(:, 2:end) +  reshape( -d, size(z(:, 2:end)));

    elseif length(  imgSize ) == 3

        d = x( 1: numel(z(1:end-1, :, :)) );
        z(1:end-1, :, :) = z(1:end-1, :, :) +reshape( d, size(z(1:end-1, :, :)));
        z(2:end, :)     = z(2:end, :) - reshape( d, size(z(2:end, :)));
        
        d = x( Np + 1 : Np + numel(z(:, 1:end-1, :)) );
        z(:, 1:end-1, :)   = z(:, 1:end-1, :) + reshape( d, size(z(:, 1:end-1, :) ));
        z(:, 2:end, :)     = z(:, 2:end, :) - reshape( d, size(z(:, 2:end, :)));
        
        d = x( Np * 2 + 1: Np * 2 + numel(z(:, :, 1:end-1)) );
        z(:, :, 1:end-1)   = z(:, :, 1:end-1) +  reshape( d, size(z(:, :, 1:end-1)));
        z(:, :, 2:end)     = z(:, :, 2:end) - reshape( d, size(z(:, :, 2:end)));

    end
    
elseif mode == 3
    
    Np = prod( imgSize );
    
    if length(  imgSize ) == 2
        
        z = ones( [ Np * 2, 1], 'single' );
        
    elseif length(  imgSize ) == 3
        
        z = ones( [ Np * 3, 1], 'single' );
        
    end
end


end