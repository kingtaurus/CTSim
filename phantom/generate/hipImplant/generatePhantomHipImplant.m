close all;

%% Load image and padd size


fprintf('Loading image ... \n');

I = imread('hip_image.jpg');
J = rgb2gray(I);
J = J';

n = [1024 704];
padsize = ( n - size(J) )/2;
J = padarray(J, padsize);

% Smoothing
fprintf('Smoothing... \n');
h = fspecial('gaussian', 10, 5);
J = medfilt2( J, [7 7]);
J = imfilter(J, h, 'replicate');


% Resize
fprintf('Resizing... \n');
factors = [3/2 3/2];
meta.DimSize = size(J) .* factors;
meta.ElementSpacing = [0.3 0.3] ./ factors;
I = resize(J, meta.DimSize);


% Segmentation
fprintf('Segmentation... '); tic;

threshAir = 50;
threshAdiposeTissue = 105;
threshSoftTissue = 200;
threshBone = 252;

maskAir = logical(I < threshAir);
maskAdiposeTissue = logical(I >= threshAir & I < threshAdiposeTissue);
maskSoftTissue = logical(I >= threshAdiposeTissue & I < threshSoftTissue );
maskBone = logical(I >= threshSoftTissue & I < threshBone);
maskCorticalBone = logical(I >= threshBone);

se = strel('disk',2);
%maskSoftTissue = imerode( maskSoftTissue, se );
maskSoftTissue = imdilate( maskSoftTissue, se );

I = zeros( size(I));
I(maskAdiposeTissue) = 1;
I(maskSoftTissue)   = 2;
I(maskBone)         = 3;
I(maskCorticalBone) = 4;

%% Remove Small Features

% Add Implant
fprintf('Painting fillings... \n'); 

ellipsesParams = [ 
    -12     -60     16      16      0       5;
    -24     -75     12      12      0       2;     
     0      -75     12      12      0       2;
    -12     -60     12      12      0       2;
    -12     -80     12      12      0       2;
    -12     -62     9       9       0       5;
     ];

J = addEllipses(I, meta, ellipsesParams);
I = max( I, J);

figure;
imagesc(I'), colormap gray;

%% Add Tissue Patterns

fprintf('Painting tissue patterns... \n');
ellipsesParams = [
     45     5   20   18     0  2;
     -15    5   18   20     0  2;
	% line pattern between fillings
	   -8     5   1   11     10  1;
	  -11.5   5   1   11     10  1;
	  -15     5   1   11     10  1;
	  -18.5   5   1   11     10  1;
	  -22     5   1   11     10  1;
    % line pattern away from fillings
	  38     8   1   10    -42  1;
	  41.5  6.5  1   10    -42  1;
	  45     5   1   10    -42  1;
	  48.5  3.5   1  10    -42  1;
	  52     2   1   10    -42  1;
      
];

I = addEllipses(I, meta, ellipsesParams);

figure;
imagesc(I'), colormap gray;
writeMetaImage(uint8(I), 'HipImplantsRealistic.mhd', meta);

%%

J = zeros( [size(I,1) size(I,2) 24], 'single' );

meta.DimSize = size(J);
meta.ElementSpacing = [meta.ElementSpacing 2];

for i = 1:24
   J(:,:,i) = I; 
end

writeMetaImage(uint8(J), 'HipImplantsRealistic-3d.mhd', meta);



return;
%% Also Save Low-Res Version

fprintf('Downsampling... \n');

factorXY = 1/2;
meta.DimSize = round(meta.DimSize * factorXY);
meta.ElementSpacing = meta.ElementSpacing / factorXY;

J = zeros(meta.DimSize, 'single');
J = imresize(I(:, :), factorXY, 'bilinear');

writeMetaImage(uint8(round(J)), 'HipImplantsLow.mhd', meta);


figure;
imagesc(J'), colormap gray;




