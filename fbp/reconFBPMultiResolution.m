function img = reconFBPMultiResolution( sino, geom, window, nrsl, scale, dispFilter )
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

sharpness = 8;

dim = ndims( sino );

nu = geom.detSize(1);
nv = geom.detSize(end);
du = geom.detSpacing(1);
ou = geom.detOffset(1);
nx = geom.reconSize(2);
ny = geom.reconSize(1);
dx = geom.reconSpacing(1);
dy = geom.reconSpacing(2);


noViews = geom.noViews;
SDD = geom.SDD;
SAD = geom.SAD;
betaMin = min(geom.betas);

% cosine weighting
sino = cosineWeighting( sino, geom );

u = (( -(nu-1)/2:(nu-1)/2)+ ou ) * du  ;
% ramp filer
if geom.flatPanel
    H       = du*designFilter2(du, window, nu, 1);
    gamma   = atan( - u / SDD );
    Tau     = abs( geom.betas(end) - geom.betas(1) ) - pi;
else
    H       = du*designEquiangularFilter2(du, SDD, window, nu, 1);
    gamma   = - u / SDD ;
    Tau     = abs( geom.betas(end) - geom.betas(1) ) - pi;
end

dbeta = abs(geom.betas(end) - geom.betas(1)) / (noViews-1);

% compute the local Nyquist freq
x = ( -(nx-1)/2:(nx-1)/2) * dx;
y = ( -(ny-1)/2:(ny-1)/2) * dy;

[xx, yy] = meshgrid(x , y);
nr = ( sqrt( xx.^2 + yy.^2 ) + 10 * dx ) * SDD / SAD / du;
freqSupport = scale ./ ( nr * dbeta );


% multi-resolution filter length
totalFilterLength = length(H);
segmentsLengths = zeros( 1, nrsl );
segmentsFreq = zeros( 1, nrsl );
currentSegmentLength = 16;
for irsl = 1 : nrsl - 1
    
    segmentsLengths( irsl ) = currentSegmentLength;
    if ( nrsl - irsl ) * currentSegmentLength < ( totalFilterLength / 3 - sum( segmentsLengths ) ) && irsl > 1
        currentSegmentLength = currentSegmentLength * 2;
    end
    segmentsFreq( irsl ) = ( 2 * sum( segmentsLengths( 1 : irsl - 1 ) ) ) / totalFilterLength;
    
end
segmentsLengths( nrsl ) = totalFilterLength / 2 - sum( segmentsLengths );
segmentsFreq( nrsl ) = 2 * sum( segmentsLengths( 1 : nrsl - 1 ) ) / totalFilterLength;
segmentsFreq( 1 ) = segmentsFreq( 2 ) / 4;

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
    
    
    fr = ( ( 1 : length(H) ) / length(H) - 0.5 ) * (1/du);
    Hshow = Hmr;
    for irsl = 1:nrsl
        Hshow(:,irsl) = fftshift( Hmr(:,irsl) );
    end
    figure;
    
    subplot(221)
    plot( fr, Hshow ); axis tight; xlabel 'frequency (mm^{-1})'
    title 'Frequency split ramp filter';
        
    freqSupport =  scale / ( ( SDD / SAD / du ) * dbeta );
    H0 = zeros( size(H) );
    for irsl = 1:nrsl
        freqWeighting = 1 ./ ( 1 + exp( sharpness * ( segmentsFreq(irsl) ./ freqSupport - 1  ) ) ) ;
        H0 = H0 + freqWeighting * fftshift( Hmr(:,irsl) );
    end
    
    
    subplot(222)
    plot(fr, H0 );
    axis tight; xlabel 'frequency (mm^{-1})'
    title 'Ramp filter at the center of rotation';
    
    
    freqSupport =  scale / ( ( 50 * SDD / SAD / du ) * dbeta );
    H5 = zeros( size(H) );
    for irsl = 1:nrsl
        freqWeighting = 1 ./ ( 1 + exp( sharpness * ( segmentsFreq(irsl) ./ freqSupport - 1  ) ) ) ;
        H5 = H5 + freqWeighting * fftshift( Hmr(:,irsl) );
    end
    
    subplot(223)
    plot(fr, H5 );
    axis tight; xlabel 'frequency (mm^{-1})'
    title 'Ramp filter at 5 cm from the center of rotation';
    
    freqSupport =  scale / ( ( 100 * SDD / SAD / du ) * dbeta );
    H10 = zeros( size(H) );
    for irsl = 1:nrsl
        freqWeighting = 1 ./ ( 1 + exp( sharpness * ( segmentsFreq(irsl) ./ freqSupport - 1  ) ) ) ;
        H10 = H10 + freqWeighting * fftshift( Hmr(:,irsl) );
    end
    
    subplot(224)
    plot(fr, H10 );
    axis tight; xlabel 'frequency (mm^{-1})'
    title 'Ramp filter at 10 cm from the center of rotation';
    
    img = 0;
    return;
end


if ~geom.shortScan
    dbeta = dbeta / 2;
end

sinoFilt = sino;
img    = zeros( geom.reconSize, 'single');

for irsl = 1 : nrsl
    
    H = Hmr(:,irsl);
    
    for view = 1:noViews
        
        % get current view's sinogram
        if dim == 3
            proj = sino(:, :, view);
        else
            proj = sino(:, view)';
        end
        
        if geom.shortScan
            %weightShortScan =  shortScanSliverWeight( gamma, geom.betas(view) - betaMin,  Tau);
            weightShortScan = shortScanParkerWeight( gamma, geom.betas(view) - betaMin, Tau );
            %weightShortScan = shortScanModifiedParkerWeight( gamma, geom.betas(view) - betaMin, Tau );
            
            if dim == 3
                weightShortScan = repmat( weightShortScan, nv, 1 );
            end
            proj = proj .* weightShortScan;
            
        end
        
        % ramp filter
        proj = filterFreq( proj, H ) ;
        
        if dim == 3
            sinoFilt(:, :,  view) = proj;
        else
            sinoFilt(:, view) = proj;
        end
        
    end
    
    
    %Scale reconstructed attenuation image to (1/cm)
    imgFreq = backProjectMex(sinoFilt, geom,  1, 0, 'back,hpd' ) * dbeta * 10;
    
    % compute local merge wieght of backprojection in different frequencies
    freqWeighting = 1 ./ ( 1 + exp( sharpness * ( segmentsFreq(irsl) ./ freqSupport - 1  ) ) ) ;
    %freqWeighting( irsl / nrsl < freqSupport  ) = 1;
    freqWeighting( ~geom.map ) = 0;
    %imdisp( freqWeighting )
    
    for iz = 1:size( img, 3)
        img(:,:,iz) = img(:,:,iz) + freqWeighting .*  imgFreq(:,:,iz);
    end
    
    
end


end
