function y = huberIsoPenalty( x, mode, delta )
% compute huber penalty function cost value, derivative, and curvature
%   input:
%       x       - image to be penalized 2D or 3D in ([nx ny], [nx ny nz])
%       mode    - ( 0 cost function )
%                 ( 1 derivative )
%                 ( 2 curvature )
%       delta   - parameter for huber function (default 0.01)
%  output:
%       y       - results
%
% Meng Wu at Stanford University
% 2012 - 2013

if nargin < 3
    delta = 0.001;
end

dim = length( size(x));

if dim == 2
    
    if mode == 0
        
        y = 0;
        y = y + sumHuber( x(1:end-1, :) -  x(2:end, :), delta );
        y = y + sumHuber( x(:, 1:end-1) - x(:, 2:end), delta );
        
    elseif mode == 1
        
        y = zeros( size(x), 'single');
        
        d = dHuber( x(1:end-1, :) -  x(2:end, :), delta);
        y( 1:end-1, :)  = y( 1:end-1, :) + d;
        y( 2:end, :)    = y( 2:end, :) - d;
        
        d = dHuber( x(:, 1:end-1) -  x(:, 2:end), delta);
        y(:, 1:end-1)   = y(:, 1:end-1) + d;
        y(:, 2:end)     = y(:, 2:end) - d;
        
        
    elseif mode == 2
        
        y = zeros( size(x), 'single');
        
        d = cHuber( x(1:end-1, :) -  x(2:end, :), delta);
        y( 1:end-1, :)  = y( 1:end-1, :) + d;
        y( 2:end, :)    = y( 2:end, :) + d;
        
        d = cHuber( x(:, 1:end-1) -  x(:, 2:end), delta);
        y(:, 1:end-1)   = y(:, 1:end-1) + d;
        y(:, 2:end)     = y(:, 2:end) + d;
        
        
    elseif mode == 3

        y = zeros( size(x), 'single');
        
        c = cHuber( x(1:end-1, :) -  x(2:end, :), delta) .* ( x(1:end-1, :) +  x(2:end, :));
        y(1:end-1, :)    = y(1:end-1, :) + c;
        y(2:end, :)      = y(2:end, :) + c;
        
        c = cHuber(  x(:, 1:end-1) -  x(:, 2:end), delta) .* ( x(:, 1:end-1) +  x(:, 2:end));
        y(:, 1:end-1)    = y(:, 1:end-1) + c;
        y(:, 2:end)      = y(:, 2:end) + c;
        
        y = y ./ x;

    else
        error('wrong mode. \n');
    end
    
elseif dim == 3
    
    if mode == 0
        
        y = 0;
        y = y + sumHuber( x(1:end-1, :, :) - x(2:end, :, :), delta );
        y = y + sumHuber( x(:, 1:end-1, :) - x(:, 2:end, :), delta );
        y = y + sumHuber( x(:, :, 1:end-1) - x(:, :, 2:end), delta );
        
    elseif mode == 1
        
        y = zeros( size(x), 'single');
        
        d = dHuber( x(1:end-1, :, :) -  x(2:end, :, :), delta);
        y(1:end-1, :, :)    = y(1:end-1, :, :) + d;
        y(2:end, :, :)      = y(2:end, :, :) - d;
        
        d = dHuber( x(:, 1:end-1, :) -  x(:, 2:end, :), delta);
        y(:, 1:end-1, :)    = y(:, 1:end-1, :) + d;
        y(:, 2:end, :)      = y(:, 2:end, :) - d;
        
        d = dHuber( x(:, :,1:end-1) - x(:, :, 2:end), delta);
        y(:, :,1:end-1)     = y(:, :,1:end-1)  + d;
        y(:, :, 2:end)      = y(:, :, 2:end)   - d;
        
    elseif mode == 2
        
        y = zeros( size(x), 'single');
        
        d = cHuber( x(1:end-1, :, :) -  x(2:end, :, :), delta);
        y(1:end-1, :, :)    = y(1:end-1, :, :) + d;
        y(2:end, :, :)      = y(2:end, :, :) + d;
        
        d = cHuber( x(:, 1:end-1, :) -  x(:, 2:end, :), delta);
        y(:, 1:end-1, :)    = y(:, 1:end-1, :) + d;
        y(:, 2:end, :)      = y(:, 2:end, :) + d;
        
        d = cHuber( x(:, :,1:end-1) - x(:, :, 2:end), delta);
        y(:, :,1:end-1)     = y(:, :,1:end-1) + d;
        y(:, :, 2:end)      = y(:, :, 2:end)  + d;
        
        
    elseif mode == 3
        
        y = zeros( size(x), 'single');

        
        c =  ( x(1:end-1, :, :) +  x(2:end, :, :));
        y(1:end-1, :, :)    = y(1:end-1, :, :) + c;
        y(2:end, :, :)      = y(2:end, :, :) + c;
        
        c = ( x(:, 1:end-1, :) +  x(:, 2:end, :));
        y(:, 1:end-1, :)    = y(:, 1:end-1, :) + c;
        y(:, 2:end, :)      = y(:, 2:end, :) + c;

        c = ( x(:, :,1:end-1) + x(:, :, 2:end));
        y(:, :,1:end-1)     = y(:, :,1:end-1) + c;
        y(:, :, 2:end)      = y(:, :, 2:end)  + c;
        
        y = y ./ x;

    else
        error('wrong penalty mode. \n');
    end
    
    
end


end

function y = sumHuber( d, delta )
d = abs( d );
temp = d.^2 / 2 ;
temp( d > delta ) = delta * d( d > delta ) - delta^2/2;
y = sum(temp(:));
end

function d = dHuber( d, delta)
d( d > delta )  = delta;
d( d < -delta ) = - delta;
end

function d = cHuber( d, delta)
qmap = ( abs(d) <= delta );
d( qmap ) = 1;
d( ~qmap ) = 0;
end
