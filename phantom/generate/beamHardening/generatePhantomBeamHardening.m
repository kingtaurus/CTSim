
n = 640;
I = zeros( n, n );
meta.DimSize = size(I) ;
meta.ElementSpacing = [0.4 0.4];

fprintf('Painting tissue patterns... \n');
ellipsesParams = [
    0       0   120  120    0   1;
    % line pattern away from fillings
   +80     0    20   20    0   2
   -80     0    20   20    0   2
     0    +80   20   20    0   2
     0    -80   20   20    0   2];

I = addEllipses(I, meta, ellipsesParams);


writeMetaImage(uint8(I), 'BeamHardening.mhd', meta);

figure;
imagesc(I'), colormap gray; axis image;


J = zeros( [n n 128] );
meta.DimSize = size(J) ;
meta.ElementSpacing = [0.4 0.4 1];
for i = 1 : size(J, 3)
   J(:,:,i) = I;
end

writeMetaImage(uint8(J), 'BeamHardening-3d.mhd', meta);


figure;
imagesc(J(:,:,end)'), colormap gray; axis image;
