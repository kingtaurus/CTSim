%% Crop

fprintf('Cropping... '); tic;

[I meta] = readMetaImage('Head+Neck_noartifacts_512x512x225_LittleEndian.mhd');
I(:, :, 1:11) = 0;
I(:, 359:end, 1:107) = 0;
I(:, 344:end, 108:125) = 0;
I(:, 316:end, 126:170) = 0;
I(:, 344:end, 171:185) = 0;
I(:, 379:end, 186:end) = 0;
% I = I(141:360, 71:360, 11:170);
I = I(141:360, 71:360, 112:156);
meta.DimSize = size(I);
writeMetaImage(I, 'Head+Neck_noartifacts_512x512x225_LittleEndian_crop.mhd', meta);
clear meta I;

fprintf('done in %0.1fs.\n', toc);


%% Extract Slices

fprintf('Extracting slices... '); tic;
[I meta] = readMetaImage('Head+Neck_noartifacts_512x512x225_LittleEndian_crop.mhd');
% lower tooth row at 43.5mm (upper tooth row would be at 28.5mm?)
x = round(21/meta.ElementSpacing(1):159/meta.ElementSpacing(1));
y = round(15/meta.ElementSpacing(2):196/meta.ElementSpacing(2));
zcenter = round(43.5/meta.ElementSpacing(3));
zwidth = ceil(20/meta.ElementSpacing(3));
z = (zcenter-zwidth):(zcenter+zwidth);
I = I(x, y, z);
meta.DimSize = size(I);

fprintf('done in %0.1fs.\n', toc);


%% Resize

fprintf('Resizing... '); tic;
factors = [4 4 1];
meta.DimSize = meta.DimSize .* factors;
meta.ElementSpacing = meta.ElementSpacing ./ factors;
I = resize(I, meta.DimSize);

fprintf('done in %0.1fs.\n', toc);


%% Segmentation

fprintf('Segmentation... '); tic;
threshAir = 600;
threshAdiposeTissue = 950;
threshSoftTissue = 1500;
threshBone = 2000;

maskAir = logical(I < threshAir);
maskAdiposeTissue = logical(I >= threshAir & I < threshAdiposeTissue);
maskSoftTissue = logical(I >= threshAdiposeTissue & I < threshSoftTissue);
maskBone = logical(I >= threshSoftTissue & I < threshBone);
maskCorticalBone = logical(I >= threshBone);

I(maskAir) = 0;
I(maskAdiposeTissue) = 1;
I(maskSoftTissue) = 2;
I(maskBone) = 3;
I(maskCorticalBone) = 4;

fprintf('done in %0.1fs.\n', toc);


%% Remove Small Features

fprintf('Filtering... '); tic;

I = medfilt3(I, 9);
writeMetaImage(uint8(I), 'Head+Neck_noartifacts_512x512x225_LittleEndian-seg.mhd', meta);

fprintf('done in %0.1fs.\n', toc);


%% Replace Outline Shape

fprintf('Outline correction... '); tic;

bigEllipseCenter = ([size(I, 1), size(I, 2)]+1) / 2;
bigEllipseAxes = [62, 84]./meta.ElementSpacing(1:2);
smallEllipseCenter = [(size(I, 1)+1)/2, 0.45*(size(I, 2)+1)];
smallEllipseAxes = [50, 68]./meta.ElementSpacing(1:2);
for x = 1:size(I, 1)
	for y = 1:size(I, 2)
		if ~isInsideEllipse([x y], bigEllipseCenter, bigEllipseAxes)
			I(x, y, :) = 0; % set outside to air
		elseif ~isInsideEllipse([x, y], smallEllipseCenter, smallEllipseAxes)
			for z = 1:size(I, 3)
				I(x, y, z) = max(I(x, y, z), 1); % fill inside with tissue (or keep whatever other matter is already there)
			end
		end
	end
end

fprintf('done in %0.1fs.\n', toc);


%% Add Fillings

fprintf('Painting fillings... '); tic;

origin = -(meta.DimSize-1)/2 .* meta.ElementSpacing;

% large filling, middle part
bigFillingMCenter = ([27.8, -39.2, 1.9] - origin)./meta.ElementSpacing;
bigFillingMRadius = [3.5, 1.2, 1.5]./meta.ElementSpacing;
bigFillingMAngle = 70/180*pi;
I = paintEllipse(I, bigFillingMCenter, bigFillingMRadius, bigFillingMAngle, 5);

% large filling, front part
bigFillingFCenter = ([26.9, -40.9, 1.4] - origin)./meta.ElementSpacing;
bigFillingFRadius = [3, 1.2, 2]./meta.ElementSpacing;
bigFillingFAngle = -20/180*pi;
I = paintEllipse(I, bigFillingFCenter, bigFillingFRadius, bigFillingFAngle, 5);

% small filling
smallFillingCenter = ([-17.6, -51.7, 1.9] - origin)./meta.ElementSpacing;
smallFillingRadius = [1.5, 1.0, 1.5]./meta.ElementSpacing;
%smallFillingRadius = [2.0, 1.5, 2.0]./meta.ElementSpacing;
smallFillingAngle = 17/180*pi;
I = paintEllipse(I, smallFillingCenter, smallFillingRadius, smallFillingAngle, 6);

fprintf('done in %0.1fs.\n', toc);


writeMetaImage(I, 'Head+Neck_noartifacts_512x512x225_fillings.mhd', meta);

%% Add Tissue Patterns

[I, meta] = readMetaImage('Head+Neck_noartifacts_512x512x225_fillings.mhd');


fprintf('Painting tissue patterns... '); tic;

origin = -(meta.DimSize-1)/2 .* meta.ElementSpacing;

lineRadius = [7, 0.65, 9]./meta.ElementSpacing;

connectingLine = bigFillingMCenter-smallFillingCenter;
lineAngle = atan2(connectingLine(2), connectingLine(1));

patternCenter(1:2) = [0 -46.5]./ meta.ElementSpacing(1:2) + (meta.DimSize(1:2)-1)/2 ;
patternCenter(3) = ceil(meta.DimSize(3)/2) * 1.6 / meta.ElementSpacing(3);

for l = -2:2
	lineCenter = patternCenter + [0, l * 2.4 /meta.ElementSpacing(2), 0];
	I = paintEllipse(I, lineCenter, lineRadius, lineAngle, 1);
end

patternCenter(1:2) = [0 36.5]./ meta.ElementSpacing(1:2)+ (meta.DimSize(1:2)-1)/2 ;
for l = -2:2
	lineCenter = patternCenter + [l * 2.4 /meta.ElementSpacing(2), 0,  0];
	I = paintEllipse(I, lineCenter, lineRadius, -1.25 , 1);
end

writeMetaImage(I, 'Head+Neck_noartifacts_512x512x225_LittleEndian-patt.mhd', meta);

fprintf('done in %0.1fs.\n', toc);


%% Save volume for CT 3D Testing

[I, meta] = readMetaImage('Head+Neck_noartifacts_512x512x225_LittleEndian-patt.mhd');

fprintf('Saving Slice for CT Testing... '); tic;
writeMetaImage(I, 'headRealistic2-3d.mhd', meta);

%% Also Save Low-Res Version

fprintf('Downsampling... '); tic;
I = single(I);
factorXY = 1/2;
meta.DimSize(1:2) = round(meta.DimSize(1:2) * factorXY);
meta.ElementSpacing(1:2) = meta.ElementSpacing(1:2) / factorXY;
J = zeros(meta.DimSize, 'single');
for z = 1:size(J, 3)
	J(:, :, z) = imresize(I(:, :, z), factorXY, 'bilinear');
end

J = uint8(round(J));
meta.DimSize =  size(J);
writeMetaImage(J, 'headRealisticLow-3d.mhd', meta);

fprintf('done in %0.1fs.\n', toc);


%% Save Single Slice for CT Testing
[I, meta] = readMetaImage('Head+Neck_noartifacts_512x512x225_LittleEndian-patt.mhd');

s = ceil(meta.DimSize(3)/2);
meta.NDims = 2;
meta.DimSize = meta.DimSize(1:2);
meta.ElementSpacing = meta.ElementSpacing(1:2);
writeMetaImage(I(:, :, s), 'headRealistic.mhd', meta);

fprintf('done in %0.1fs.\n', toc);


%% Save Top and Bottom Slices for Dose Simulations

fprintf('Saving Top and Bottom Slices for Dose Simulations... '); tic;

[I, meta] = readMetaImage('Head+Neck_noartifacts_512x512x225_LittleEndian-patt.mhd');

s = ceil(meta.DimSize(3)/2);
meta.NDims = 2;
meta.DimSize = meta.DimSize(1:2);
meta.ElementSpacing = meta.ElementSpacing(1:2);
writeMetaImage(I(:, :, s-3), 'headRealistic-TopSlice.mhd', meta);

meta.NDims = 2;
meta.DimSize = meta.DimSize(1:2);
meta.ElementSpacing = meta.ElementSpacing(1:2);
writeMetaImage(I(:, :, s+3), 'headRealistic-BottomSlice.mhd', meta);


fprintf('done in %0.1fs.\n', toc);
