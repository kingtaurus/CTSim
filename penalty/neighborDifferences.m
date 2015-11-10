function y = neighborDifferences( x, mode )

if nargin < 2
    mode = 0;
end

dim = length( size(x));

axy = 1 / sqrt(2);
az = 1 / 2;
y = zeros( size(x), 'single');

if dim == 2
    
    if mode == 0
        
        d =  x(1:end-1, :) -  x(2:end, :);
        y( 1:end-1, :)  = y( 1:end-1, :) + d;
        y( 2:end, :)    = y( 2:end, :) - d;
        
        d = x(:, 1:end-1) -  x(:, 2:end);
        y(:, 1:end-1)   = y(:, 1:end-1) + d;
        y(:, 2:end)     = y(:, 2:end) - d;
        
        d = x(1:end-1, 1:end-1) -  x(2:end, 2:end);
        y( 1:end-1, 1:end-1)    = y( 1:end-1, 1:end-1) + axy * d;
        y( 2:end, 2:end)        = y( 2:end, 2:end) - axy * d;
        
        d = x(2:end, 1:end-1) -  x( 1:end-1, 2:end);
        y( 2:end, 1:end-1)      = y(2:end, 1:end-1)  + axy * d;
        y( 1:end-1, 2:end)      = y( 1:end-1, 2:end) - axy * d;

    elseif mode == 1
        
        d =  x(1:end-1, :) +  x(2:end, :);
        y( 1:end-1, :)  = y( 1:end-1, :) + d;
        y( 2:end, :)    = y( 2:end, :) + d;
        
        d = x(:, 1:end-1) +  x(:, 2:end);
        y(:, 1:end-1)   = y(:, 1:end-1) + d;
        y(:, 2:end)     = y(:, 2:end) + d;
        
        d = x(1:end-1, 1:end-1) +  x(2:end, 2:end);
        y( 1:end-1, 1:end-1)    = y( 1:end-1, 1:end-1) + axy * d;
        y( 2:end, 2:end)        = y( 2:end, 2:end) + axy * d;
        
        d = x(2:end, 1:end-1) +  x( 1:end-1, 2:end);
        y( 2:end, 1:end-1)      = y(2:end, 1:end-1)  + axy * d;
        y( 1:end-1, 2:end)      = y( 1:end-1, 2:end) + axy * d;
        
        
    else
        error('wrong mode. \n');
    end
    
elseif dim == 3
    
    if mode == 0

        d = x(1:end-1, :, :) -  x(2:end, :, :);
        y(1:end-1, :, :)    = y(1:end-1, :, :) + d;
        y(2:end, :, :)      = y(2:end, :, :) - d;
        
        d = x(:, 1:end-1, :) -  x(:, 2:end, :);
        y(:, 1:end-1, :)    = y(:, 1:end-1, :) + d;
        y(:, 2:end, :)      = y(:, 2:end, :) - d;
        
        d = x(1:end-1, 1:end-1, :) -  x(2:end, 2:end, :);
        y(1:end-1, 1:end-1, :)  = y(1:end-1, 1:end-1, :) + axy * d;
        y(2:end, 2:end, :)      = y(2:end, 2:end, :) - axy * d;
        
        d = x(2:end, 1:end-1, :) -  x(1:end-1, 2:end, :);
        y(2:end, 1:end-1, :)    = y(2:end, 1:end-1, :) + axy * d;
        y(1:end-1, 2:end, :)    = y(1:end-1, 2:end, :) - axy * d;
        
        d = x(:, :,1:end-1) - x(:, :, 2:end);
        y(:, :,1:end-1)     = y(:, :,1:end-1)  + az * d;
        y(:, :, 2:end)      = y(:, :, 2:end)   - az * d;

    elseif mode == 1
        
        d = x(1:end-1, :, :) +  x(2:end, :, :);
        y(1:end-1, :, :)    = y(1:end-1, :, :) + d;
        y(2:end, :, :)      = y(2:end, :, :) + d;
        
        d = x(:, 1:end-1, :) -  x(:, 2:end, :);
        y(:, 1:end-1, :)    = y(:, 1:end-1, :) + d;
        y(:, 2:end, :)      = y(:, 2:end, :) + d;
        
        d = x(1:end-1, 1:end-1, :) -  x(2:end, 2:end, :);
        y(1:end-1, 1:end-1, :)  = y(1:end-1, 1:end-1, :) + axy * d;
        y(2:end, 2:end, :)      = y(2:end, 2:end, :) + axy * d;
        
        d = x(2:end, 1:end-1, :) -  x(1:end-1, 2:end, :);
        y(2:end, 1:end-1, :)    = y(2:end, 1:end-1, :) + axy * d;
        y(1:end-1, 2:end, :)    = y(1:end-1, 2:end, :) + axy * d;
        
        d = x(:, :,1:end-1) - x(:, :, 2:end);
        y(:, :,1:end-1)     = y(:, :,1:end-1)  + az * d;
        y(:, :, 2:end)      = y(:, :, 2:end)   + az * d;
        
    else
        error('wrong mode. \n');
    end
    
else
    error('wrong dimension!\n');
end