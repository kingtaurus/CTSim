% Generate numerical bead geomerty module of the CAT phantom
%
% Meng Wu at Stanford University
% 2014

ellipsesParams = zeros(1, 6 );
ellipsesParams(1,:) = [0    0   90  90  0   1];
ellipsesParamsDensity = [0    0   90  90  0   0.95; ...
    0    0   75  75  0   1];

I = zeros( [1024 1024], 'single' );
meta.DimSize = size(I) ;
meta.ElementSpacing = [0.25 0.25];
I = addEllipses(I, meta, ellipsesParams);

J = zeros( [1024 1024], 'single' );
J = addEllipses(J, meta, ellipsesParamsDensity);

figure;
imagesc(I .* J), colormap gray; axis image;

%%
K = zeros([ 1024 1024 160], 'single' );
meta.DimSize = size(K) ;
meta.ElementSpacing = [0.25 0.25 0.25];
density = K;

a = round( 1 / meta.ElementSpacing(3) ); 
b = round( 0.25 / meta.ElementSpacing(3) );
for i = 1:meta.DimSize(3)
    
    j = round( i - meta.DimSize(3)/2 );
    S = I;
    
    % course ramps
    if mod( j, a ) == 0 &&  abs( floor( j / a ) ) < 20
        S( j * a/2 + 512, 312 ) = 100;
        S( -j * a/2 + 512, 712 ) = 100;
        
        S( 312, j * a/2 + 512 ) = 100;
        S( 712, -j * a/2 + 512 ) = 100;
        
        
    end
    
    % fine ramps
    if mod(j, b) == 0 && abs( floor( j / b )  ) <= 12
        S( j * 8 * b + 512, 362 ) = 100;
        S( -j * 8 * b + 512, 662 ) = 100;
    end
    S( 752, 512) = 100;
    
    if j == 0
        S( 362, 512) = 100;
        S( 662, 512) = 100;
    end
    
    K(:,:,i) = S;
    density(:,:,i) = J;
end


figure;
imagesc( sum(K,3), [150 170] ), colormap gray; axis image;

figure;
imagesc( squeeze( sum(K( :, 300:800,: ),2) ) ), colormap gray; axis image;

%% save to raw data
writeMetaImage(uint8(K), 'Catphan-geometry-3d.mhd', meta);
writeMetaImage(density, 'Catphan-geometry-3d-density.mhd', meta);

