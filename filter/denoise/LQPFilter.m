%% LQPFilter.m

%% LQPFilter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% uhat = LQPFilter(v,w,s)
% Apply a local quadratic programming filter with window half-width w,
% neighborhood window half-width wn, averaging weight standard
% deviation s and similarity norm standard deviation a to the image v
% Inputs:
%    v    : the input image
%    w    : the window half-width
%    s    : the standard deviation of noise
% Outputs:
%    uhat : the output image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [uhat, t] = LQPFilter(v,w,s)
% get the size of the image
options = optimset('quadprog');
options = optimset(options,'Algorithm', 'interior-point-convex',...
    'Display','off');

n = length(v);
u = BilateralFilter( v, 3, 1, s*3 );
% iterate over the pixels
ut = zeros(n);
h = waitbar(0,'message') ;
for x = 1:n
    for y = 1:n
        % compute the search window
        iMin = max(x-w,1);
        iMax = min(x+w,n);
        jMin = max(y-w,1);
        jMax = min(y+w,n);
        Ns   = v(iMin:iMax , jMin:jMax);
        
        Ms = u(iMin:iMax , jMin:jMax) - u(x,y);
        
        m = length(Ns(:));
        z = Ms(:);
        
        
        w = exp( - z.^2 / ( 2 * s^2 ) );
        z = z - 0.8 * sum( z.* w ) / sum(w);
        %z = z .* max(abs(z) - s/2, 0) ./ ( abs(z) + 1);
        
        H = z*z' + s^2 * eye(m);
        w = quadprog(H,[],[],[],ones(1,m),1,[],[],[],options );
        
        w = w / sum(w);
        
        t = t + sum(w(:).^2 ) * 20^2 + sum(  w.* Ms(:) )^2;
        
        % compute the weighted average over the search neighborhood
        ut(x,y) = sum(w(:) .* Ns(:));
        
    end
    waitbar(x / n, h);
    
end


close(h);
uhat = ut;

end