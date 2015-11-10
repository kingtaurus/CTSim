
n = 512;
I = zeros( n, n );
meta.DimSize = size(I) ;
meta.ElementSpacing = [1 1];

fprintf('Painting tissue patterns... \n');
ellipsesParams = [
     0    +55   130   130    0   1
     0    -55   130   130    0   1];

I = addEllipses(I, meta, ellipsesParams);

I( end/2-54:end/2+55, end/2-129:end/2+130 ) = 1;


figure;
imagesc(I'), colormap gray; axis image;
writeMetaImage(uint8(I), 'JangAEC.mhd', meta);

J = zeros( [n n 256] );
meta.DimSize = size(J) ;
meta.ElementSpacing = [1 1 1];
for i = 1 : size(J, 3)
   J(:,:,i) = I;
end

writeMetaImage(uint8(J), 'JangAEC-3d.mhd', meta);

figure;
imagesc(J(:,:,end)'), colormap gray; axis image;
