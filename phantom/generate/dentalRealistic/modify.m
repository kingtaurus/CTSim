[J meta] = readMetaImage('headRealistic.mhd');
origin = -(meta.DimSize-1)/2 .* meta.ElementSpacing;

figure;
imagesc(I'), colormap gray; axis image;

factors = [2 2];
meta.DimSize = size(J) .* factors;
meta.ElementSpacing = meta.ElementSpacing ./ factors;
I = imresize(J, meta.DimSize,'nearest');



fprintf('Painting tissue patterns... \n');
ellipsesParams = [
    35      0   13   13     0  2;
   -45      3   12   12     0  2;
    % line pattern between fillings
    -40     3   0.6   9    -15  1;
    -42.5   3   0.6   9    -15  1;
    -45     3   0.6   9    -15  1;
    -47.5   3   0.6   9    -15  1;
    -50     3   0.6   9    -15  1;
    % line pattern away from fillings
    30     0   0.6   9     -15  1;
    32.5   0   0.6   9     -15  1;
    35     0   0.6   9     -15  1;
    37.5   0   0.6   9     -15  1;
    40     0   0.6   9     -15  1;
    
    ];

I = addEllipses(I, meta, ellipsesParams);


figure;
imagesc(I'), colormap gray; axis image;


writeMetaImage(uint8(I), 'headRealistic.mhd', meta);

