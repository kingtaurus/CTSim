function y = totalVariationIsoPenalty( x, mode, epsilon )
% compute total variation penalty function cost value, derivative, and curvature
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


epsilon = 1e-8;

dim = length( size(x));

if dim == 2
    
    d1 = x(1:end-1, 1:end-1) - x(2:end, 1:end-1);
    d2 = x(1:end-1, 1:end-1) - x(1:end-1, 2:end);
    c  = 1 ./ sqrt( d1.^2 + d2.^2 + epsilon );
    
    if mode == 0        
        y = sum( 1 ./ c(:) );

    elseif mode == 1
        
        y = zeros( size(x), 'single');
        
        y( 1:end-1,  1:end-1)  = y( 1:end-1,  1:end-1) + d1 .* c;
        y( 2:end,  1:end-1)    = y( 2:end,  1:end-1) - d1 .* c;
        
        y(1:end-1, 1:end-1)   = y(1:end-1, 1:end-1) + d2 .* c ;
        y(1:end-1, 2:end)     = y(1:end-1, 2:end) - d2 .* c;
        
    elseif mode == 2 || mode == 3
        
        y = zeros( size(x), 'single');
        
        y( 1:end-1,  1:end-1)  = y( 1:end-1,  1:end-1) +  c;
        y( 2:end,  1:end-1)    = y( 2:end,  1:end-1) + c;
        
        y(1:end-1, 1:end-1)   = y(1:end-1, 1:end-1) + c ;
        y(1:end-1, 2:end)     = y(1:end-1, 2:end) + c;
        
    else
        error('wrong mode. \n');
    end
    
elseif dim == 3
    
    d1 = x(1:end-1, 1:end-1, 1:end-1) - x(2:end, 1:end-1, 1:end-1);
    d2 = x(1:end-1, 1:end-1, 1:end-1) - x(1:end-1, 2:end, 1:end-1);
    d3 = x(1:end-1, 1:end-1, 1:end-1) - x(1:end-1, 1:end-1, 2:end);
    c  = 1 ./ sqrt( d1.^2 + d2.^2 + d3.^2 + epsilon );

    if mode == 0
        
        y = sum( 1 ./ c(:) );
        
    elseif mode == 1
        
        y = zeros( size(x), 'single');
        
        y(1:end-1, 1:end-1, 1:end-1)    = y(1:end-1, 1:end-1, 1:end-1) + d1 .* c;
        y(2:end, 1:end-1, 1:end-1)      = y(2:end, 1:end-1, 1:end-1)   - d1 .* c;
        
        y(1:end-1, 1:end-1, 1:end-1)    = y(1:end-1, 1:end-1, 1:end-1) + d2 .* c;
        y(1:end-1, 2:end, 1:end-1)      = y(1:end-1, 2:end, 1:end-1)   - d2 .* c;
        
        y(1:end-1, 1:end-1, 1:end-1)    = y(1:end-1, 1:end-1,1:end-1)  + d3 .* c;
        y(1:end-1, 1:end-1, 2:end)      = y(1:end-1, 1:end-1, 2:end)   - d3 .* c;
        
    elseif mode == 2 || mode == 3
        
        y = zeros( size(x), 'single');
        
        y(1:end-1, 1:end-1, 1:end-1)    = y(1:end-1, 1:end-1, 1:end-1) + c;
        y(2:end, 1:end-1, 1:end-1)      = y(2:end, 1:end-1, 1:end-1)   + c;
        
        y(1:end-1, 1:end-1, 1:end-1)    = y(1:end-1, 1:end-1, 1:end-1) + c;
        y(1:end-1, 2:end, 1:end-1)      = y(1:end-1, 2:end, 1:end-1)   + c;
        
        y(1:end-1, 1:end-1, 1:end-1)    = y(1:end-1, 1:end-1,1:end-1)  + c;
        y(1:end-1, 1:end-1, 2:end)      = y(1:end-1, 1:end-1, 2:end)   + c;
        
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
