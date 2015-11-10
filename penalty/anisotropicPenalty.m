function y = anisotropicPenalty( x, mode, delta )
% compute anisotropix quadratic penalty function cost value, derivative, and curvature
%   input:
%       x       - image to be penalized 2D or 3D in ([nx ny], [nx ny nz])
%       mode    - ( 0 cost function )
%                 ( 1 derivative )
%                 ( 2 curvature )
%       delta   - parameter for huber function (default 0.1)
%  output:
%       y       - results
%
% Meng Wu at Stanford University
% 2012 - 2013

if nargin < 3
    delta = 0.1;
end

dim = length( size(x));

axy = 1 / sqrt(2);
az = 1 / 2;

if dim == 2
    
    if mode == 0
        
        y = 0;
        y = y + sumAniso( x(1:end-1, :) -  x(2:end, :), delta );
        y = y + sumAniso( x(:, 1:end-1) - x(:, 2:end), delta );
        
        y = y + axy * sumAniso( x(1:end-1, 1:end-1) - x(2:end, 2:end), delta );
        y = y + axy * sumAniso( x(2:end, 1:end-1) - x(1:end-1, 2:end), delta );
        
    elseif mode == 1
        
        y = zeros( size(x), 'single');
        
        d = dAniso( x(1:end-1, :) -  x(2:end, :), delta);
        y( 1:end-1, :)  = y( 1:end-1, :) + d;
        y( 2:end, :)    = y( 2:end, :) - d;
        
        d = dAniso( x(:, 1:end-1) -  x(:, 2:end), delta);
        y(:, 1:end-1)   = y(:, 1:end-1) + d;
        y(:, 2:end)     = y(:, 2:end) - d;
        
        d = dAniso( x(1:end-1, 1:end-1) -  x(2:end, 2:end), delta);
        y( 1:end-1, 1:end-1)    = y( 1:end-1, 1:end-1) + axy * d;
        y( 2:end, 2:end)        = y( 2:end, 2:end) - axy * d;
        
        d = dAniso( x(2:end, 1:end-1) -  x( 1:end-1, 2:end), delta);
        y( 2:end, 1:end-1)      = y(2:end, 1:end-1)  + axy * d;
        y( 1:end-1, 2:end)      = y( 1:end-1, 2:end) - axy * d;
        
    elseif mode == 2
        
        y = zeros( size(x), 'single');
        
        d = cAniso( x(1:end-1, :) -  x(2:end, :), delta);
        y( 1:end-1, :)  = y( 1:end-1, :) + d;
        y( 2:end, :)    = y( 2:end, :) + d;
        
        d = cAniso( x(:, 1:end-1) -  x(:, 2:end), delta);
        y(:, 1:end-1)   = y(:, 1:end-1) + d;
        y(:, 2:end)     = y(:, 2:end) + d;
        
        d = cAniso( x(1:end-1, 1:end-1) -  x(2:end, 2:end), delta);
        y( 1:end-1, 1:end-1)    = y( 1:end-1, 1:end-1) + axy * d;
        y( 2:end, 2:end)        = y( 2:end, 2:end) + axy * d;
        
        d = cAniso( x(2:end, 1:end-1) -  x( 1:end-1, 2:end), delta);
        y( 2:end, 1:end-1)      = y( 2:end, 1:end-1)  + axy * d;
        y( 1:end-1, 2:end)      = y( 1:end-1, 2:end) + axy * d;
        
    else
        error('wrong penalty mode. \n');
    end
    
elseif dim == 3
    
    if mode == 0
        
        y = 0;
        y = y + sumAniso( x(1:end-1, :, :) - x(2:end, :, :), delta );
        y = y + sumAniso( x(:, 1:end-1, :) - x(:, 2:end, :), delta );
        
        y = y + axy * sumAniso( x(1:end-1, 1:end-1, :) - x(2:end, 2:end, :), delta );
        y = y + axy * sumAniso( x(2:end, 1:end-1, :) - x(1:end-1, 2:end, :), delta );
        
        y = y + az * sumAniso( x(:, :, 1:end-1) - x(:, :, 2:end), delta );
        
    elseif mode == 1
        
        y = zeros( size(x), 'single');
        
        d = dAniso( x(1:end-1, :, :) -  x(2:end, :, :), delta);
        y(1:end-1, :, :)    = y(1:end-1, :, :) + d;
        y(2:end, :, :)      = y(2:end, :, :) - d;
        
        d = dAniso( x(:, 1:end-1, :) -  x(:, 2:end, :), delta);
        y(:, 1:end-1, :)    = y(:, 1:end-1, :) + d;
        y(:, 2:end, :)      = y(:, 2:end, :) - d;
        
        d = dAniso( x(1:end-1, 1:end-1, :) -  x(2:end, 2:end, :), delta);
        y(1:end-1, 1:end-1, :)  = y(1:end-1, 1:end-1, :) + axy * d;
        y(2:end, 2:end, :)      = y(2:end, 2:end, :) - axy * d;
        
        d = dAniso( x(2:end, 1:end-1, :) -  x(1:end-1, 2:end, :), delta);
        y(2:end, 1:end-1, :)    = y(2:end, 1:end-1, :)  + axy * d;
        y(1:end-1, 2:end, :)    = y(1:end-1, 2:end, :) - axy * d;
        
        d = dAniso( x(:, :,1:end-1) - x(:, :, 2:end), delta);
        y(:, :,1:end-1)     = y(:, :,1:end-1)  + az * d;
        y(:, :, 2:end)      = y(:, :, 2:end)   - az * d;
        
        
    elseif mode == 2
        
        y = zeros( size(x), 'single');
        
        d = cAniso( x(1:end-1, :, :) -  x(2:end, :, :), delta);
        y(1:end-1, :, :)    = y(1:end-1, :, :) + d;
        y(2:end, :, :)      = y(2:end, :, :) + d;
        
        d = cAniso( x(:, 1:end-1, :) -  x(:, 2:end, :), delta);
        y(:, 1:end-1, :)    = y(:, 1:end-1, :) + d;
        y(:, 2:end, :)      = y(:, 2:end, :) + d;
        
        d = cAniso( x(1:end-1, 1:end-1, :) -  x(2:end, 2:end, :), delta);
        y(1:end-1, 1:end-1, :)  = y(1:end-1, 1:end-1, :) + axy * d;
        y(2:end, 2:end, :)      = y(2:end, 2:end, :) + axy * d;
        
        d = cAniso( x(2:end, 1:end-1, :) -  x(1:end-1, 2:end, :), delta);
        y(2:end, 1:end-1, :)    = y(2:end, 1:end-1, :)  + axy * d;
        y(1:end-1, 2:end, :)    = y(1:end-1, 2:end, :) + axy * d;
        
        d = cAniso( x(:, :,1:end-1) - x(:, :, 2:end), delta);
        y(:, :,1:end-1)     = y(:, :,1:end-1)  + az * d;
        y(:, :, 2:end)      = y(:, :, 2:end)   + az * d;
       
        
    else
        error('wrong penalty mode. \n');
    end
    
end
end

function y = sumAniso( d, delta )
    d = d.^2 .* cAniso(d, delta) / 2 ;
    y = sum( d(:));
end

function d = dAniso( d, delta )
    d = d.*cAniso( d, delta );
end

function d = cAniso( d, delta )
    d = exp( - d.^2 / delta^2 );
end