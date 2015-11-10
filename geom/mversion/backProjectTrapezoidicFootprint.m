function img = backProjectTrapezoidicFootprint (sino, geom, isFlatDetector)
% Forward projection using distance-driven method in 2D
% input:
%       img - image to project
%       phan - image parameters including spacing and size
%       geom - system geometry
% output:
%       sinp - projection result aka sinogram
%
% Meng Wu
% 2011.12

if nargin < 3
    isFlatDetector = true;
end

nx = geom.reconSize(1);
ny = geom.reconSize(2);
dx = geom.reconSpacing(1);
dy = dx;

du = geom.detSpacing;
nu = geom.detSize;

img = zeros(nx, ny);

wx = - (nx-1)/2.;
wy = - (ny-1)/2.;
wu = - (nu-1)/2.;

dsd = geom.SDD;
dso = geom.SAD;
noViews = geom.noViews;
dsd_ds = dsd / du;

for view = 1:noViews
    
    %compute the projection angle
    beta = -geom.betas(view);
    
    %compute one project view
    oneView = sino(:, view) * dx / 10;
    
    sin_a = sin(beta);
    cos_a = cos(beta);
    
    dx_sin_2 = dx * sin_a / 2.;
    dx_cos_2 = dx * cos_a / 2.;
    dy_sin_2 = dy * sin_a / 2.;
    n_dy_cos_2 = -dy * cos_a / 2.;
    
    
    %this can be futher accelerate, but complicated ToT
    
    for iy = 1 : ny
        for ix = 1 : nx

            
            xc = (ix - 1 + wx) * dx; % center
            yc = (iy - 1 + wy) * dy;
            
            tproj0 = xc * cos_a - yc * sin_a;
            dsxy0 = dso + xc * sin_a + yc * cos_a;
            tproj0_dsxy0 = sqrt(tproj0 * tproj0 + dsxy0 * dsxy0);
            
            lu = amplitudeFunction(dsxy0, tproj0, tproj0_dsxy0, sin_a, cos_a);
            
            if isFlatDetector
                taus = boundariesTrapezoidicFootPrintFlat(dsd_ds, tproj0, dsxy0, wu, dx_cos_2, dx_sin_2, dy_sin_2, n_dy_cos_2);
            else
                taus = boundariesTrapezoidicFootPrintArc(dsd_ds, tproj0, dsxy0, wu, dx_cos_2, dx_sin_2, dy_sin_2, n_dy_cos_2);
            end
            
            [ weight_u, isproj, iumin, iumax ] = weightTrapezoidicFootPrint(taus, nu);
            
            %if outside the detection region
            if ~isproj
                continue;
            end
            
            temp = 0;
            iiu = 1;
            for iu = iumin : iumax
                %add to projection view
                temp = temp + oneView(iu+1) * weight_u(iiu);
                iiu = iiu + 1;
            end
            
            img(ix, iy) = img(ix, iy) + temp * lu;
            
        end %iy
        
    end %ix
    
    
    
end

end






