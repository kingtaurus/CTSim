
ellipsesParams = [0    0   100  65  0   1;
    0   0   50  50   0   2];
ellipsesParamsDensity = [0    0   100  65  0  1;
                         0   0   50  50   0   0.75];

I = zeros( [512 512], 'single' );
meta2.DimSize = size(I) ;
meta2.ElementSpacing = [0.5 0.5];
I1 = addEllipses(I, meta2, ellipsesParams);
I2 = addEllipses(I, meta2, ellipsesParams(1,:));

J = zeros( [512 512], 'single' );
J1 = addEllipses(J, meta2, ellipsesParamsDensity);
J2 = addEllipses(J, meta2, ellipsesParamsDensity(1,:));

figure;
imagesc(I1 .* J1, [0 2]), colormap gray; axis image;

figure;
imagesc(I2 .* J2, [0 2]), colormap gray; axis image;

%%
K = zeros([ 512 512 128], 'single' );
meta.DimSize = size(K) ;
meta.ElementSpacing = [0.5 0.5 0.5];
density = K;

for i = 1:meta.DimSize(3)
    
    if mod( i - 10, 20 ) >= 10
        K(:,:,i) = I1;     
        density(:,:,i) = J1;
    else
        K(:,:,i) = I2;
        density(:,:,i) = J2;
    end
    
    
    
    imagesc(K(:,:,i), [0 2]), colormap gray; axis image;
    
end

%%
figure;
imagesc( K(:,:,end/2)); colormap gray; axis image;

return;

%% save to raw data
writeMetaImage(uint8(K), 'HBP2-phantom-3d.mhd', meta);
writeMetaImage(density, 'HBP2-phantom-3d-density.mhd', meta);

