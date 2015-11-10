function y = quadPenalty( x, mode )
% compute quadratic penalty function cost value, derivative, and curvature
%   input:
%       x       - image to be penalized 2D or 3D in ([nx ny], [nx ny nz])
%       mode    - ( 0 cost function )
%                 ( 1 derivative )
%                 ( 2 curvature )
%  output:
%       y       - results
%
% Meng Wu at Stanford University
% 2012 - 2013


dim = length( size(x));
axy = 1 / sqrt(2);
az = 1 / 2;
    
if dim == 2

    if mode == 0
        
        y = 0;
        y = y + sumQuad( x(1:end-1, :) -  x(2:end, :));
        y = y + sumQuad( x(:, 1:end-1) - x(:, 2:end));
        
        y = y + axy * sumQuad( x(1:end-1, 1:end-1) - x(2:end, 2:end));
        y = y + axy * sumQuad( x(2:end, 1:end-1) - x(1:end-1, 2:end));
        
    elseif mode == 1
        
        y = zeros( size(x), 'single');
        
        d = dQuad( x(1:end-1, :) -  x(2:end, :));
        y( 1:end-1, :)  = y( 1:end-1, :) + d;
        y( 2:end, :)    = y( 2:end, :) - d;
        
        d = dQuad( x(:, 1:end-1) -  x(:, 2:end));
        y(:, 1:end-1)   = y(:, 1:end-1) + d;
        y(:, 2:end)     = y(:, 2:end) - d;
        
        d = dQuad( x(1:end-1, 1:end-1) -  x(2:end, 2:end));
        y( 1:end-1, 1:end-1)    = y( 1:end-1, 1:end-1) + axy * d;
        y( 2:end, 2:end)        = y( 2:end, 2:end) - axy * d;
        
        d = dQuad( x(2:end, 1:end-1) -  x( 1:end-1, 2:end));
        y( 2:end, 1:end-1)      = y( 2:end, 1:end-1)  + axy * d;
        y( 1:end-1, 2:end)      = y( 1:end-1, 2:end) - axy * d;
        
    elseif mode == 2
        
    y = ones( size(x), 'single') * ( 4 + 4*axy);
        
    else
        error('wrong penalty mode. \n');
    end
    
elseif dim == 3
     if mode == 0
        
        y = 0;
        y = y + sumQuad( x(1:end-1, :, :) - x(2:end, :, :));
        y = y + sumQuad( x(:, 1:end-1, :) - x(:, 2:end, :));
        
        y = y + axy * sumQuad( x(1:end-1, 1:end-1, :) - x(2:end, 2:end, :));
        y = y + axy * sumQuad( x(2:end, 1:end-1, :) - x(1:end-1, 2:end, :));
        
        y = y + az * sumQuad( x(:, :, 1:end-1) - x(:, :, 2:end));

    elseif mode == 1
        
        y = zeros( size(x), 'single');
        
        d = dQuad( x(1:end-1, :, :) -  x(2:end, :, :));
        y(1:end-1, :, :)    = y(1:end-1, :, :) + d;
        y(2:end, :, :)      = y(2:end, :, :) - d;
        
        d = dQuad( x(:, 1:end-1, :) -  x(:, 2:end, :));
        y(:, 1:end-1, :)    = y(:, 1:end-1, :) + d;
        y(:, 2:end, :)      = y(:, 2:end, :) - d;
        
        d = dQuad( x(1:end-1, 1:end-1, :) -  x(2:end, 2:end, :));
        y(1:end-1, 1:end-1, :)  = y(1:end-1, 1:end-1, :) + axy * d;
        y(2:end, 2:end, :)      = y(2:end, 2:end, :) - axy * d;
        
        d = dQuad( x(2:end, 1:end-1, :) -  x(1:end-1, 2:end, :));
        y(2:end, 1:end-1, :)    = y(2:end, 1:end-1, :)  + axy * d;
        y(1:end-1, 2:end, :)    = y(1:end-1, 2:end, :) - axy * d;
        
        d = dQuad( x(:, :,1:end-1) - x(:, :, 2:end));
        y(:, :,1:end-1)     = y(:, :,1:end-1)  + az * d;
        y(:, :, 2:end)      = y(:, :, 2:end)   - az * d;
        
        
    elseif mode == 2
        
        y = ones( size(x), 'single') * ( 4 + 4*axy + 2*az);
        
    else
        error('wrong penalty mode. \n');
    end

end

end

function y = sumQuad( d )
    d = d.^2 / 2 ;
    y = sum( d(:));
end

function d = dQuad( d )

end

