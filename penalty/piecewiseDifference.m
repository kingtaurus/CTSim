function z = piecewiseDifference( x, imgSize,  mode )




if mode == 1
    
    Np = numel( x );
    
    if length(  imgSize ) == 2
        
        z = zeros( [ Np * 4, 1], 'single' );
        
        d = x(1:end-1, :) -  x(2:end, :);
        z( 1: numel(d) ) = d(:);
        
        d = x(:, 1:end-1) -  x(:, 2:end);
        z( Np + 1: Np + numel(d) ) = d(:);
        
        d = x(1:end-1, 1:end-1) -  x(2:end, 2:end);
        z( Np * 2 + 1: Np * 2 + numel(d) ) = d(:);
        
        d = x(2:end, 1:end-1) -  x( 1:end-1, 2:end);
        z( Np * 3 + 1: Np * 3 + numel(d) ) = d(:);

    elseif length(  imgSize ) == 3
        
        z = zeros( [Np * 5, 1], 'single' );
        
        d = x(1:end-1, :, :) - x(2:end, :, :);
        z( 1: numel(d) ) = d(:);
        
        d =  x(:, 1:end-1, :) - x(:, 2:end, :);
        z( Np + 1: Np + numel(d) ) = d(:);
        
        d = x(1:end-1, 1:end-1, :) - x(2:end, 2:end, :);
        z( Np * 2 + 1: Np * 2 + numel(d) ) = d(:);
        
        d = x(2:end, 1:end-1, :) - x(1:end-1, 2:end, :);
        z( Np * 3 + 1: Np * 3 + numel(d) ) = d(:);
        
        d = x(:, :, 1:end-1) - x(:, :, 2:end);
        z( Np * 4 + 1: Np * 4 + numel(d) ) = d(:);
        
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
        
        d = x( Np * 2 + 1: Np * 2 +  numel(z(1:end-1, 1:end-1)) );
        z(1:end-1, 1:end-1) = z(1:end-1, 1:end-1) + reshape( d, size(z(1:end-1, 1:end-1)));
        z(2:end, 2:end)     = z(2:end, 2:end) + reshape( -d, size(z(2:end, 2:end)));
        
        d = x( Np * 3 + 1: Np * 3 + numel(z(2:end, 1:end-1)) );
        z(2:end, 1:end-1)   = z(2:end, 1:end-1)+ reshape( d, size(z(2:end, 1:end-1) ));
        z( 1:end-1, 2:end)  = z( 1:end-1, 2:end) + reshape( -d, size(z( 1:end-1, 2:end) ));
        
        
    elseif length(  imgSize ) == 3
        
        
        d = x( 1: numel(z(1:end-1, :, :)) );
        z(1:end-1, :, :) = z(1:end-1, :, :) +reshape( d, size(z(1:end-1, :, :)));
        z(2:end, :)     = z(2:end, :) - reshape( d, size(z(2:end, :)));
        
        d = x( Np + 1 : Np + numel(z(:, 1:end-1, :)) );
        z(:, 1:end-1, :)   = z(:, 1:end-1, :) + reshape( d, size(z(:, 1:end-1, :) ));
        z(:, 2:end, :)     = z(:, 2:end, :) - reshape( d, size(z(:, 2:end, :)));
        
        d = x( Np * 2 + 1: Np * 2 +  numel(z(1:end-1, 1:end-1, :)) );
        z(1:end-1, 1:end-1, :) = z(1:end-1, 1:end-1, :) + reshape( d, size(z(1:end-1, 1:end-1, :)));
        z(2:end, 2:end, :)     = z(2:end, 2:end, :) - reshape( d, size(z(2:end, 2:end, :)));
        
        d = x( Np * 3 + 1: Np * 3 + numel(z(2:end, 1:end-1, :)) );
        z(2:end, 1:end-1, :)   = z(2:end, 1:end-1, :) + reshape( d, size(z(2:end, 1:end-1, :)));
        z(1:end-1, 2:end, :)   = z(1:end-1, 2:end, :) - reshape( d, size(z(1:end-1, 2:end, :)));
        
        d = x( Np * 4 + 1: Np * 4 + numel(z(:, :, 1:end-1)) );
        z(:, :, 1:end-1)   = z(:, :, 1:end-1) +  reshape( d, size(z(:, :, 1:end-1)));
        z(:, :, 2:end)     = z(:, :, 2:end) - reshape( d, size(z(:, :, 2:end)));
        
        
    end
    
elseif mode == 3
    
    axy = 1 / sqrt(2);
    az = 1 / 2;
    
    Np = prod( imgSize );
    
    if length(  imgSize ) == 2
        z = zeros( [ Np * 4, 1], 'single' );
        z( 1: Np * 2 ) = 1;
        z( Np * 2 + 1: Np * 4 ) = axy;
        
    elseif length(  imgSize ) == 3
        
        z = zeros( [ Np * 5, 1], 'single' );
        z( 1: Np * 2 ) = 1;
        z( Np * 2 + 1: Np * 4 ) = axy;
        z( Np * 4 + 1: Np * 5 ) = az;
        
        
    end
end


end