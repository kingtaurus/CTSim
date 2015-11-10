function p = ellipses(nx, ny, ellipsesParams, dx, dy)
% generate 2D image with ellipeses parameters
% intputs:
%       nx, ny
%       ellipsesParams
%       dx, dy
% outputs:
%       p


p = zeros(ny,nx);
xx = [(-nx+1)/2 : (nx-1)/2] * dx;
yy = [(-ny+1)/2 : (ny-1)/2] * dy;

[ xg, yg ] = meshgrid( yy, xx ); 


for k = 1:size(ellipsesParams,1)
    
    x0 = ellipsesParams(k, 1);           % x offset
    y0 = ellipsesParams(k, 2);           % y offset
    
    asq = ellipsesParams(k, 3)^2;        % a^2
    bsq = ellipsesParams(k, 4)^2;        % b^2
    
    phi = ellipsesParams(k, 5)*pi/180;   % first Euler angle in radians
    
    A = ellipsesParams(k, 6);            % Amplitude change for this ellipsoid
    
    x=xg-x0;                             % Center the ellipsesParams
    y= yg -y0;
    cosp = cos(phi);
    sinp = sin(phi);
    
    t = ((x.*cosp + y.*sinp).^2)./asq + ((y.*cosp - x.*sinp).^2)./bsq <= 1 ;
    p(t) = p(t) + A;
    
end



end