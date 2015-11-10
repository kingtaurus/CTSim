% Generate numerical sensitivity module of the CAT phantom
%
% Meng Wu at Stanford University
% 2014

supraRadius = 58.4;
supraTargetsRaduis = 6.1;
supraAngle = [pi/12 pi/3 pi/2 2*pi/3 11*pi/12];


subRadius = 15;
subDensity = 1.1;
subTargetsDiameters = [10 8 6 4 3 2 1.5 1 0];
subTargetsRaduis = subTargetsDiameters / 2;
subDistancs = 7;

ellipsesParams = zeros(12, 6 );
ellipsesParams(1,:) = [0    0   90  90  0   1];
ellipsesParamsDensity = [0    0   90  90  0   0.95; ...
                        0    0   75  75  0   1];

t = 3;
for i = 1:2
    if i == 1
        alpha = 0;
    else
        alpha = pi;
    end
    
    for j = 1:5
        ellipsesParams((i-1)*5+j+1,:) =  [ supraRadius*cos( alpha + supraAngle(j) ) ...
           -supraRadius*sin( alpha + supraAngle(j) )  supraTargetsRaduis  supraTargetsRaduis    0   t];
        t = t + 1;
    end
    
end


alpha = pi / 12 ;
for j = 1:8
        ellipsesParams((i-1)*8+j+12,:) =  [ subRadius*cos( alpha ) ...
            subRadius*sin( alpha )  subTargetsRaduis(j)  subTargetsRaduis(j)    0   t+1];
    alpha = alpha + ( subTargetsRaduis(j) + subTargetsRaduis(j+1) + subDistancs ) / subRadius;
end


I = zeros( [1024 1024], 'single' );
meta.DimSize = size(I) ;
meta.ElementSpacing = [0.2 0.2];
I = addEllipses(I, meta, ellipsesParams);


J = zeros( [1024 1024], 'single' );
J = addEllipses(J, meta, ellipsesParamsDensity);


figure;
imagesc(I .* J), colormap gray; axis image;

return;

%%


K = zeros([ 1024 1024 100], 'single' );
meta.DimSize = size(K) ;
meta.ElementSpacing = [0.2 0.2 0.4];

density = K;
for i = 1:meta.DimSize(3)
    K(:,:,i) = I;
    density(:,:,i) = J;
end

%% save to raw data
writeMetaImage(uint8(K), 'Catphan-sensitivity-3d.mhd', meta);
writeMetaImage(density, 'Catphan-sensitivity-3d-density.mhd', meta);

