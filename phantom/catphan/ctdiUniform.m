ellipsesParams = zeros(2, 6 );
ellipsesParams(1,:) = [0   0   160  160    0   1];

I = zeros( [1024 1024], 'single' );
meta.DimSize = size(I) ;
meta.ElementSpacing = [0.4 0.4];
I = addEllipses(I, meta, ellipsesParams );

K = I;
for i = 1 : 44
    K = K + imrotate( I, i , 'crop' ); 
end
I = K / 45;


J = ones([ 1024 1024 160], 'single' );
meta.DimSize = size(J) ;
meta.ElementSpacing = [0.4 0.4 0.5];

density = J;
for i = 1:meta.DimSize(3)
    %J(:,:,i) = single( I > 0 );
    density(:,:,i) = I;
end

figure; imdisp( I );

%%
meta.DimSize = size(I) ;
meta.ElementSpacing = [0.4 0.4];
writeMetaImage(uint8(I), 'CTDI-uniform.mhd', meta);
writeMetaImage(I, 'CTDI-uniform-density.mhd', meta);


%% save to raw data
meta.DimSize = size(J) ;
meta.ElementSpacing = [0.4 0.4 0.5];
writeMetaImage(uint8(J), 'CTDI-uniform-3d.mhd', meta);
writeMetaImage(density, 'CTDI-uniform-3d-density.mhd', meta);

