I = zeros( 256,  256);
meta.DimSize = size(I) ;
meta.ElementSpacing = [1.2 1.2];

fprintf('Painting tissue patterns... \n');
ellipsesParams = [
    0       0    120   100   0   3;
    0       0    115    95   0   1
    % line pattern away from fillings
    %     +30     +30   15   15    0   2
    %     -30     +30   15   15    0   2
    %     +30     -30   15   15    0   2
    -0     -0   90   80    0   2
    -0     -0   90   80    0   4    ];

I = addEllipses(I, meta, ellipsesParams);

figure;
imagesc(I'), colormap gray; axis image;


writeMetaImage(uint8(I), 'dynamicHead2.mhd', meta);