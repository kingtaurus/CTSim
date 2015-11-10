ellipsesParams = zeros(2, 6 );
ellipsesParams(1,:) = [0   0   90  90    0   1];

I = zeros( [1024 1024], 'single' );
meta.DimSize = size(I) ;
meta.ElementSpacing = [0.2 0.2];
I = addEllipses(I, meta, ellipsesParams );

J = zeros([ 1024 1024 160], 'single' );
meta.DimSize = size(J) ;
meta.ElementSpacing = [0.2 0.2 0.25];

density = J;
for i = 1:meta.DimSize(3)
    J(:,:,i) = single( I > 0 );
    density(:,:,i) = I;
end

figure; imdisp( I );

%%
meta.DimSize = size(I) ;
meta.ElementSpacing = [0.2 0.2];
writeMetaImage(uint8(I), 'Catphan-uniform.mhd', meta);
writeMetaImage(I, 'Catphan-uniform-density.mhd', meta);


%% save to raw data
meta.DimSize = size(J) ;
meta.ElementSpacing = [0.2 0.2 0.25];
writeMetaImage(uint8(J), 'Catphan-uniform-3d.mhd', meta);
writeMetaImage(density, 'Catphan-uniform-3d-density.mhd', meta);

