function s = measurePathImages( x, t, type, d )


if type == 1
    s = std( x(:) - t(:) ) * 1000 / 0.202 ;
elseif type == 2
    s = mean( abs( x(:) - t(:) ) ) * 1000 / 0.202 ;
elseif type == 3
    s = std( x(:) - t(:) ) / std( d(:)) ;
elseif type == 4
    s = mean( abs( x(:) - t(:) ) ) / mean( abs( d(:)) );
elseif type == 5
    
    x = x( d(:) > 0.001 );
    t = t( d(:) > 0.001 );
    d = d( d(:) > 0.001 );
    s = mean( abs( x(:) - t(:) ) ./ d(:) );
end