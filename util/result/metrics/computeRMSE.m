function e = computeRMSE(img1, img2, roi)

if nargin < 2
	error('Two images have to be specified!');
end
if ndims(img1) ~= ndims(img2)
	error('The given images do not have an equal number of dimensions!');
end
if ~all(size(img1) == size(img2))
	error('The given images do not have the same size!');
end
[M N] = size(img1);

if nargin < 3
   roi = ones(M, N);
end
if isvector(roi) && length(roi) == 4
	roiextents = roi;
	roi = zeros(M, N);
	roi(roiextents(1):roiextents(2), roiextents(3):roiextents(4)) = 1;
end
if (ndims(roi) ~= 2) || ~all(size(roi) == [M N])
	error('roi has to be either a 4-vector of rectangular extents or a mask!');
end
roi = logical(roi);

% compute root-mean-square error
diff = img1-img2;
diff = diff(roi);
diff = diff(:);
diff = diff - mean(diff);
e = sqrt(mean(diff.^2));

%fprintf('The rmse is %3i HU. \n', round(e) );

