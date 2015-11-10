function y = hyperbolaPenalty( x, mode, delta )
% compute hyperbola penalty function cost value, derivative, and curvature
%   input:
%       x       - image to be penalized 2D or 3D in ([nx ny], [nx ny nz])
%       mode    - ( 0 cost function )
%                 ( 1 derivative )
%                 ( 2 curvature )
%       delta   - parameter for hyperbola function (default 0.001)
%  output:
%       y       - results
%
% Meng Wu at Stanford University
% 2012 - 2013

if nargin < 3
    delta = 0.001;
end

dim = length( size(x));
axy = 1 / sqrt(2);
az = 1 / 2;
    
if dim == 2

    if mode == 0
        
        y = 0;
        y = y + sumHyper( x(1:end-1, :) -  x(2:end, :), delta );
        y = y + sumHyper( x(:, 1:end-1) - x(:, 2:end), delta );
        
        y = y + axy * sumHyper( x(1:end-1, 1:end-1) - x(2:end, 2:end), delta );
        y = y + axy * sumHyper( x(2:end, 1:end-1) - x(1:end-1, 2:end), delta );
        
    elseif mode == 1
        
        y = zeros( size(x), 'single');
        
        d = dHyper( x(1:end-1, :) -  x(2:end, :), delta);
        y( 1:end-1, :)  = y( 1:end-1, :) + d;
        y( 2:end, :)    = y( 2:end, :) - d;
        
        d = dHyper( x(:, 1:end-1) -  x(:, 2:end), delta);
        y(:, 1:end-1)   = y(:, 1:end-1) + d;
        y(:, 2:end)     = y(:, 2:end) - d;
        
        d = dHyper( x(1:end-1, 1:end-1) -  x(2:end, 2:end), delta);
        y( 1:end-1, 1:end-1)    = y( 1:end-1, 1:end-1) + axy * d;
        y( 2:end, 2:end)        = y( 2:end, 2:end) - axy * d;
        
        d = dHyper( x(2:end, 1:end-1) -  x( 1:end-1, 2:end), delta);
        y( 2:end, 1:end-1)      = y(2:end, 1:end-1)  + axy * d;
        y( 1:end-1, 2:end)      = y( 1:end-1, 2:end) - axy * d;
        
    elseif mode == 2
        
        y = zeros( size(x), 'single');
        
        d = cHyper( x(1:end-1, :) -  x(2:end, :), delta);
        y( 1:end-1, :)  = y( 1:end-1, :) + d;
        y( 2:end, :)    = y( 2:end, :) + d;
        
        d = cHyper( x(:, 1:end-1) -  x(:, 2:end), delta);
        y(:, 1:end-1)   = y(:, 1:end-1) + d;
        y(:, 2:end)     = y(:, 2:end) + d;
        
        d = cHyper( x(1:end-1, 1:end-1) -  x(2:end, 2:end), delta);
        y( 1:end-1, 1:end-1)    = y( 1:end-1, 1:end-1) + axy * d;
        y( 2:end, 2:end)        = y( 2:end, 2:end) + axy * d;
        
        d = cHyper( x(2:end, 1:end-1) -  x( 1:end-1, 2:end), delta);
        y( 2:end, 1:end-1)      = y(2:end, 1:end-1)  + axy * d;
        y( 1:end-1, 2:end)      = y( 1:end-1, 2:end) + axy * d;
        
    else
        error('wrong penalty mode. \n');
    end
    
elseif dim == 3
     if mode == 0
        
        y = 0;
        y = y + sumHyper( x(1:end-1, :, :) - x(2:end, :, :), delta );
        y = y + sumHyper( x(:, 1:end-1, :) - x(:, 2:end, :), delta );
        
        y = y + axy * sumHyper( x(1:end-1, 1:end-1, :) - x(2:end, 2:end, :), delta );
        y = y + axy * sumHyper( x(2:end, 1:end-1, :) - x(1:end-1, 2:end, :), delta );
        
        y = y + az * sumHyper( x(:, :, 1:end-1) - x(:, :, 2:end), delta );

    elseif mode == 1
        
        y = zeros( size(x), 'single');
        
        d = dHyper( x(1:end-1, :, :) -  x(2:end, :, :), delta);
        y(1:end-1, :, :)    = y(1:end-1, :, :) + d;
        y(2:end, :, :)      = y(2:end, :, :) - d;
        
        d = dHyper( x(:, 1:end-1, :) -  x(:, 2:end, :), delta);
        y(:, 1:end-1, :)    = y(:, 1:end-1, :) + d;
        y(:, 2:end, :)      = y(:, 2:end, :) - d;
        
        d = dHyper( x(1:end-1, 1:end-1, :) -  x(2:end, 2:end, :), delta);
        y(1:end-1, 1:end-1, :)  = y(1:end-1, 1:end-1, :) + axy * d;
        y(2:end, 2:end, :)      = y(2:end, 2:end, :) - axy * d;
        
        d = dHyper( x(2:end, 1:end-1, :) -  x(1:end-1, 2:end, :), delta);
        y(2:end, 1:end-1, :)    = y(2:end, 1:end-1, :)  + axy * d;
        y(1:end-1, 2:end, :)    = y(1:end-1, 2:end, :) - axy * d;
        
        d = dHyper( x(:, :,1:end-1) - x(:, :, 2:end), delta);
        y(:, :,1:end-1)     = y(:, :,1:end-1)  + az * d;
        y(:, :, 2:end)      = y(:, :, 2:end)    - az * d;
        
    elseif mode == 2
        
        y = zeros( size(x), 'single');
        
        d = cHyper( x(1:end-1, :, :) -  x(2:end, :, :), delta);
        y(1:end-1, :, :)    = y(1:end-1, :, :) + d;
        y(2:end, :, :)      = y(2:end, :, :) + d;
        
        d = cHyper( x(:, 1:end-1, :) -  x(:, 2:end, :), delta);
        y(:, 1:end-1, :)    = y(:, 1:end-1, :) + d;
        y(:, 2:end, :)      = y(:, 2:end, :) + d;
        
        d = cHyper( x(1:end-1, 1:end-1, :) -  x(2:end, 2:end, :), delta);
        y(1:end-1, 1:end-1, :)  = y(1:end-1, 1:end-1, :) + axy * d;
        y(2:end, 2:end, :)      = y(2:end, 2:end, :) + axy * d;
        
        d = cHyper( x(2:end, 1:end-1, :) -  x(1:end-1, 2:end, :), delta);
        y(2:end, 1:end-1, :)    = y(2:end, 1:end-1, :)  + axy * d;
        y(1:end-1, 2:end, :)    = y(1:end-1, 2:end, :) + axy * d;
        
        d = cHyper( x(:, :,1:end-1) - x(:, :, 2:end), delta);
        y(:, :,1:end-1)     = y(:, :,1:end-1)  + az * d;
        y(:, :, 2:end)      = y(:, :, 2:end)   + az * d;
        
    else
        error('wrong penalty mode. \n');
    end

end

end


function y = sumHyper( d, delta )
    y = delta^2 * sum( sqrt( 1 + (d(:)/delta).^2 ) - 1 );
end

function d = dHyper( d, delta)
    d = d ./ sqrt( 1 + (d/delta).^2 );
end

function d = cHyper( d, delta)
    d = 1 ./ ( sqrt( 1 + (d/delta).^2 ).^3 );
end
