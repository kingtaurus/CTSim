function sino = forwardProjectRayDriven( img, phan, geom )
% Forward projection using Ray-driven method in 2D
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

overSamplingDepth = 1;

du = geom.detSpacing;
nu = geom.detSize;
dx = phan.spacing(1);

SAD = geom.SAD;
ADD = geom.ADD;
noViews = geom.noViews;

sino = zeros(nu, noViews);

cubemin = ( -(phan.size-1)/2 + phan.offset ) .* dx;
cubemax = (  (phan.size-1)/2 + phan.offset ) .* dx;
wxy     = - (phan.size-1)/2. + phan.offset ;
wu      = - (nu-1)/2. + geom.detOffset(1);

% forward-project views
for view = 1:noViews
    
    % camera angle between the negative y axis and center to camera position in clockwise direction
    beta = geom.betas(view);
    
    % determine sampling distance from minimal pixel spacing
    deltaSampling = min(dx)/overSamplingDepth;
    
    % source position
    sourcePos = -SAD * [sin(beta), cos(beta)];
    
    % init this view's sinogram
    sinoView = NaN(nu, 1);
    
    % cast a ray for each detector pixel
    for d = 1:nu
        
        % angle from principle ray to current ray on detector
        gamma = atan( (d - 1 + wu ) * du / (SAD+ADD) );
        
        % combined angle between the negative y axis and current ray
        phi = beta + gamma;
        
        % ray direction
        rayDir = [sin(phi), cos(phi)];
        
        % accumulate attenuation values
        proj = 0;
        [hit, tn, tf] = intersectRayWithCuboid(sourcePos, rayDir, cubemin, cubemax);
        if hit
            
            sampleDepths = tn*(1+10*eps):deltaSampling:tf/(1+10*eps);
            noDepths = length(sampleDepths);
            sampleIndexes = zeros(noDepths, 2);
            for i = 1:noDepths
                % calculate sampling point
                samplePos = sourcePos + sampleDepths(i)*rayDir;
                
                % interpolate sample
                sampleIndexes(i, :) = samplePos./dx - wxy + 1; % Matlab indexes are 1-based!
            end
            % interpolate sample using Matlab's built-in interp2 method
            proj = sum(interp2(img, sampleIndexes(:, 2), sampleIndexes(:, 1), '*linear'));
        end
        
        % multiply with sampling step size to obtain projection integrals
        sinoView(d) = proj*(deltaSampling/10); % convert deltaSampling to cm before multiplying it with the attenuation sum that is in 1/cm
    end
    
    % downsample intermediate projections that were oversampled
    sino(:, view) = sinoView;
    
end

end
