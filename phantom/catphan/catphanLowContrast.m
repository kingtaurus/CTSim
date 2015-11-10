supraRadius = 50;


supraDensitis = [ 1.1 1.05 1.03];
supraTargetsDiameters = [15 9 8 7 6 5 4 3 2 1];
supraTargetsRaduis = supraTargetsDiameters / 2;
supraDistancs = 4;

subRadius = 25;
subDensity = 1.1;
subTargetsDiameters = [9 7 5 3 1];
subTargetsRaduis = subTargetsDiameters / 2;
subDistancs = 4;

ellipsesParams = zeros(41, 6 );
ellipsesParams(1,:) = [0   0   90  90    0   0.95];
ellipsesParams(2,:) = [0   0   75  75    0   1];

alpha = 0;
for j = 1:9
    for i = 1:3
        ellipsesParams((i-1)*9+j+2,:) =  [ supraRadius*cos( alpha + i*2*pi/3 ) ...
            supraRadius*sin( alpha + i*2*pi/3 )  supraTargetsRaduis(j)  supraTargetsRaduis(j)    0   supraDensitis(i)];
    end
    alpha = alpha + ( supraTargetsRaduis(j) + supraTargetsRaduis(j+1) + supraDistancs ) / supraRadius;
end


for i = 1:3
    alpha = pi / 12 ;
    for j = 1:4
        ellipsesParams((i-1)*4+j+29,:) =  [ subRadius*cos( alpha + i*2*pi/3 ) ...
            subRadius*sin( alpha + i*2*pi/3 )  subTargetsRaduis(j)  subTargetsRaduis(j)    0   subDensity];
        alpha = alpha + ( subTargetsRaduis(j) + subTargetsRaduis(j+1) + subDistancs ) / subRadius;
    end
    
end

I = zeros( [1024 1024], 'single' );
meta.DimSize = size(I) ;
meta.ElementSpacing = [0.2 0.2];
I = addEllipses(I, meta, ellipsesParams(1:2,:));

I0 = zeros( [1024 1024], 'single' );
I0 = addEllipses(I0, meta, ellipsesParams(1:29,:));

I1 = zeros( [1024 1024], 'single' );
I1 = addEllipses(I1, meta, ellipsesParams(1:33,:));

I2 = zeros( [1024 1024], 'single' );
I2 = addEllipses(I2, meta, ellipsesParams(1:37,:));

I3 = zeros( [1024 1024], 'single' );
I3 = addEllipses(I3, meta, ellipsesParams(1:41,:));


figure;
imagesc(I0', [0.8 1.4]), colormap gray; axis image;

figure;
imagesc(I1', [0.8 1.4]), colormap gray; axis image;

figure;
imagesc(I2', [0.8 1.4]), colormap gray; axis image;

figure;
imagesc(I3', [0.8 1.4]), colormap gray; axis image;



J = zeros([ 1024 1024 160], 'single' );
meta.DimSize = size(J) ;
meta.ElementSpacing = [0.2 0.2 0.25];

density = J;
for i = 1:meta.DimSize(3)
    J(:,:,i) = single( I > 0 );
    
    if abs( i - meta.DimSize(3)/2 + 1/2 ) * meta.ElementSpacing(3) <= 1.5
        density(:,:,i) = I3;
    elseif abs( i - meta.DimSize(3)/2 + 1/2 ) * meta.ElementSpacing(3) <= 2.5
        density(:,:,i) = I2;
    elseif abs( i - meta.DimSize(3)/2 + 1/2 ) * meta.ElementSpacing(3) <= 3.5
        density(:,:,i) = I1;
    else
        density(:,:,i) = I0;
    end
    
end

%% save to raw data
writeMetaImage(uint8(J), 'Catphan-contrast-3d.mhd', meta);
writeMetaImage(density, 'Catphan-contrast-3d-density.mhd', meta);

