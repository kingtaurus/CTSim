% Generate numerical bead geomerty module of the CAT phantom
%
% Meng Wu at Stanford University
% 2014

ellipsesParams = zeros(51, 6 );
ellipsesParams(1,:) = [0    0   90  90  0   1];
ellipsesParamsDensity = [0    0   90  90  0   0.95; ...
    0    0   75  75  0   1];


radius = 52;
patternLength = 10;

distance = 15;
GapSizes = [2, 1.6, 1.4, 1.2, 1.1, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3 0.2];

n = 2;
alpha = 0;
for i =  1 : length( GapSizes )
    
    gs = GapSizes(i);
    
    t = 5;
        
    if i <= 2
        t = 3;
    elseif i <= 4
        t = 4;
    end
    
    alpha0 =  ( alpha + atan(  ( gs * t ) / radius ) ) * 180 / pi;
    
    for j = 1:t
        
        ellipsesParams(n,:) =  [ radius*cos( alpha ) ...
            -radius*sin( alpha )  patternLength*2  gs/2  -alpha0   2];
        
        n = n + 1;
        alpha = alpha + 2 * atan( gs / radius );
        
        
    end
    
    alpha = alpha + 2 * atan(  distance / radius / 2);
    
    
    
end
ellipsesInert = [0 0 radius-patternLength/2 radius-patternLength/2 0 1 ];
ellipsesOutter = [0 0 radius+patternLength/2 radius+patternLength/2 0 1 ];


I = zeros( [1024 1024], 'single' );
meta.DimSize = size(I) ;
meta.ElementSpacing = [0.2 0.2];
I = addEllipses(I, meta, ellipsesParams);


N = zeros( [1024 1024], 'single' );
N = addEllipses(N, meta, ellipsesInert);
I(N == 1) = 1;

N = zeros( [1024 1024], 'single' );
N = addEllipses(N, meta, ellipsesOutter);
I( N == 0 & I > 0 ) = 1;


N = zeros( [1024 1024], 'single' );
ellipsesParamsFlat = [0    0   90  90  0   1];
N = addEllipses(N, meta, ellipsesParamsFlat);

J = zeros( [1024 1024], 'single' );
J = addEllipses(J, meta, ellipsesParamsDensity);

figure;
imagesc( I ), colormap gray; axis image;

%%
K = zeros([ 1024 1024 20], 'single' );
meta.DimSize = size(K) ;
meta.ElementSpacing = [0.2 0.2 0.5];
density = K;
for i = 1:meta.DimSize(3)
    
    j = round( i - meta.DimSize(3)/2 );
    
    if abs( j ) <= 2
        K(:,:,i) = I;
    else
        K(:,:,i) = N;
    end
    density(:,:,i) = J;
    
end
figure;
imagesc( sum(K,3) ), colormap gray; axis image;

%% save to raw data
writeMetaImage(uint8(K), 'Catphan-resolution-3d.mhd', meta);
writeMetaImage(density, 'Catphan-resolution-3d-density.mhd', meta);


%%
meta.ElementSpacing = [0.4 0.4 0.5];
writeMetaImage(uint8(K), 'Catphan-resolution-x2-3d.mhd', meta);
writeMetaImage(density, 'Catphan-resolution-x2-3d-density.mhd', meta);


