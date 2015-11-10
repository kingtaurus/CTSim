function [img, imgMultFreq] = reconFBP2dMultiResolution( sino, geom, window, nrsl, scale, dispFilter )
% Fan-beam CT filtered backpeoject reconstruction
% multi-resoltion piecewise linear frequency response FBP
% inout:
%       sino    - log sinogram
%       geom    - geometry parameters
%       window  - window function
%       nrsl    - number of resolution component
%       scale   - scale of frequency support ( default 0.5 )
% output:
%       img  - reconstructed attenuation image (1/cm)
%
% Meng Wu, Stanford University and University Erlangen-Nuremberg
% 2014

if nargin < 6
    dispFilter = false;
end


nu = geom.detSize;
du = geom.detSpacing;
noViews = geom.noViews;
SDD = geom.SDD;
SAD = geom.SAD;
minbetas = min( geom.betas );

nx = geom.reconSize(2);
ny = geom.reconSize(1);
dx = geom.reconSpacing(1);
dy = geom.reconSpacing(1);

x = ( -(nx-1)/2:(nx-1)/2) * dx;
y = ( -(ny-1)/2:(ny-1)/2) * dy;

dbeta = abs(geom.betas(end) - geom.betas(1)) / (noViews-1);
[xx, yy] = meshgrid(x , y);
nr = sqrt( xx.^2 + yy.^2 ) * SDD / SAD / du;
freqSupport = scale ./ ( nr * dbeta );


% array of detector positions in mm
u = (( -(nu-1)/2:(nu-1)/2)+ geom.detOffset ) * du  ;

% cosine weighting factors and filter
if geom.flatPanel
    weight  = SAD ./ sqrt( u.^2 + SDD^2 );
    H       = du*designFilter2(du, window, nu, 1);
    gamma   = atan( - u / SDD );
    Tau     = abs( geom.betas(end) - geom.betas(1) ) - pi;
else
    weight  = SAD / SDD * cos( u / SDD );
    H       = du*designEquiangularFilter2(du, SDD, window, nu, 1);
    gamma   = - u / SDD ;
    Tau     = abs( geom.betas(end) - geom.betas(1) ) - pi;
end

% multi-resolution filter length
totalFilterLength = length(H);
segmentsLengths = zeros( 1, nrsl );
segmentsFreq = zeros( 1, nrsl );
currentSegmentLength = 8;
for irsl = 1 : nrsl - 1
    
    segmentsLengths( irsl ) = currentSegmentLength;
    if ( nrsl - irsl ) * currentSegmentLength < ( totalFilterLength / 2 - sum( segmentsLengths ) ) && irsl > 1
        currentSegmentLength = currentSegmentLength * 2;
    end
    segmentsFreq( irsl ) = ( 2 * sum( segmentsLengths( 1 : irsl - 1 ) ) ) / totalFilterLength;
    
end
segmentsLengths( nrsl ) = totalFilterLength / 2 - sum( segmentsLengths );
segmentsFreq( nrsl ) = 2 * sum( segmentsLengths( 1 : nrsl - 1 ) ) / totalFilterLength;

% design multi-resolution filter
Hmr = zeros( length(H), nrsl );

for irsl = 1 : nrsl
    if irsl == 1
        L = segmentsLengths(1);
        Hmr( 1 : L, irsl ) = 1;
    elseif irsl == 2
        L = segmentsLengths(2);
        S = segmentsLengths(1);
        Hmr( S + 1 : S + L, irsl ) =  (L:-1:1)/L;
    else
        L1 = segmentsLengths(irsl-1);
        L2 = segmentsLengths(irsl);
        S1 = sum( segmentsLengths(1:irsl-2) );
        Hmr( S1+1 : S1 + L1 , irsl ) =  (1:L1)/L1;
        Hmr( S1 + L1 + 1 : S1 + L1 + L2 , irsl ) =  (L2:-1:1)/L2;
    end
    Hmr(:,irsl) = Hmr(:,irsl) + flipud( Hmr(:,irsl) );
    Hmr(:,irsl) = Hmr(:,irsl) .* H;
end

% Display the filters
if dispFilter
    Hshow = Hmr;
    for irsl = 1:nrsl
        Hshow(:,irsl) = fftshift( Hmr(:,irsl) );
    end
    figure;
    plot( Hshow ); axis tight;
    title 'Frequency spilt ramp filter';
end


imgMultFreq       = cell( nrsl, 1);
img    = zeros( geom.reconSize, 'single');

factor = 1;
if ~geom.shortScan
    factor = factor / 2;
end

sinoFilt = sino;
for irsl = 1 : nrsl
    
    H = Hmr(:,irsl);
    
    % filter projections
    for view = 1:noViews
        % get current view's sinogram
        proj = sino(:, view)';
        
        if geom.shortScan
            %short scan weighting
            weightShortScan = shortScanSliverWeight( gamma, geom.betas(view) - minbetas,  Tau);
            % cosine weighting
            proj = weight .* weightShortScan .* proj;
        else
            % cosine weighting
            proj = weight .* proj;
        end
        % get current view's sinogram
        sinoFilt(:, view) = filterFreq( proj, H ) ;
        
    end
    
    %Scale reconstructed attenuation image to (1/cm)
    imgFreq = backProjectMex(sinoFilt, geom,  1, 0, 'back,pd' ) * dbeta * 10 * factor;
    
    % compute local merge wieght of backprojection in different frequencies
    if irsl == 1
        freqWeighting = 1  ;
    else
        freqWeighting = 1 ./ ( 1 + exp( 6 * ( segmentsFreq(irsl) ./ freqSupport - 1  ) ) ) ;
        %freqWeighting( irsl / nrsl < freqSupport  ) = 1;
    end
    
    %imdisp( freqWeighting )
    
    img = img + freqWeighting .*  imgFreq;
    
    imgMultFreq{irsl} = imgFreq;
    
end

img( ~ geom.map ) = 0;

end


