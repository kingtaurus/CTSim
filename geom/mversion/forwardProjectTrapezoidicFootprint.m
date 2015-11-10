function sino = forwardProjectTrapezoidicFootprint(img, phan, geom, isFlatDetector)
% Forward projection using trapezoidic footprint method in 2D % Forward projection using distance-driven method in 2D
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
    geom = phan;
    phan.spacing = geom.reconSpacing;
    phan.size    = geom.reconSize;
    phan.offset  = geom.reconOffset;
end


if nargin < 4
    isFlatDetector = true;
end

du = geom.detSpacing;
nu = geom.detSize;
dx = phan.spacing(1);
dy = dx;

[nx ny] = size(img);

wx = - (nx-1)/2. + geom.reconOffset(1);
wy = - (ny-1)/2. + geom.reconOffset(2);
wu = - (nu-1)/2. + geom.detOffset(1);

dsd = geom.SDD;
dso = geom.SAD;
noViews = geom.noViews;
dsd_ds = dsd / du;

sino = zeros(nu, noViews);

for view = 1:noViews
    
    %compute the projection angle
    beta = geom.betas(view);
    
    %compute one project view
    oneView = zeros(nu, 1);
    
    sin_a = sin(beta);
    cos_a = cos(beta);
    
    dx_sin_2 = dx * sin_a / 2.;
    dx_cos_2 = dx * cos_a / 2.;
    dy_sin_2 = dy * sin_a / 2.;
    n_dy_cos_2 = -dy * cos_a / 2.;
    
    
    %this can be futher accelerate, but complicated ToT
    for ix = 1 : nx
        
        for iy = 1 : ny
            
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
            
            voxelValue = lu * img(ix, iy);
            iiu = 1;
            for iu = iumin : iumax
                %add to projection view
                oneView(iu+1) = oneView(iu+1) + voxelValue * weight_u(iiu);
                iiu = iiu + 1;
            end
            
        end %iy
        
    end %ix
    
    %compute one project view
    sino(:, view) = oneView / 10 * dx; %covert to cm
    
    
end

end




