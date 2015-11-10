function img = backProjectPixelDriven (sino, geom )
% Backprojection using pixel driven method
%
% Meng Wu
% 2011.12

du = geom.detSpacing;
nu = geom.detSize;
dx = geom.reconSpacing(1);
nx = geom.reconSize(1);
ny = geom.reconSize(2);

SAD = geom.SAD;
SDD = geom.SDD;
noViews = geom.noViews;

img = zeros(nx, ny);

wx = - (nx-1)/2. + geom.reconOffset(1) / dx;
wy = - (ny-1)/2. + geom.reconOffset(2) / dx;
wu = - (nu-1)/2. + geom.detOffset(1) / du;

% backprojection to all image pixels
for view = 1:noViews
    
    proj = sino(:, view);
    
    %compute the projection angle
    beta = geom.betas(view);
    
    sin_a = sin(beta);
    cos_a = cos(beta);
    
    for iy = 1 : ny
        
        for ix = 1 : nx
            
            xc = (ix + wx - 1) * dx; % center
            yc = (iy + wy - 1) * dx;
            xRot = xc * cos_a - yc * sin_a;
            yRot = xc * sin_a + yc * cos_a;
            
            w2 = SDD / ( SAD + yRot );
            
            u = xRot * w2 / du - wu + 1;
            
            if u >= 1 && u <= nu
                iu = floor(u);
                delta = u - iu;
                img(ix, iy) =  img(ix, iy) + w2^2 * ((1-delta)*proj(iu) + delta*proj(iu+1));
            end
            
        end %ix
        
    end %iy
    
end

end



