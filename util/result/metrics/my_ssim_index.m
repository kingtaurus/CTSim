function [mssim, ssim_map] = my_ssim_index(img1, img2, roi, K, window, L)

%========================================================================
%SSIM Index, Version 1.0
%Copyright(c) 2003 Zhou Wang
%All Rights Reserved.
%
%The author is with Howard Hughes Medical Institute, and Laboratory
%for Computational Vision at Center for Neural Science and Courant
%Institute of Mathematical Sciences, New York University.
%
%----------------------------------------------------------------------
%Permission to use, copy, or modify this software and its documentation
%for educational and research purposes only and without fee is hereby
%granted, provided that this copyright notice and the original authors'
%names appear on all copies and supporting documentation. This program
%shall not be used, rewritten, or adapted as the basis of a commercial
%software or hardware product without first obtaining permission of the
%authors. The authors make no representations about the suitability of
%this software for any purpose. It is provided "as is" without express
%or implied warranty.
%----------------------------------------------------------------------
%
%This is an implementation of the algorithm for calculating the
%Structural SIMilarity (SSIM) index between two images. Please refer
%to the following paper:
%
%Z. Wang, A. C. Bovik, H. R. Sheikh, and E. P. Simoncelli, "Image
%quality assessment: From error measurement to structural similarity"
%IEEE Transactios on Image Processing, vol. 13, no. 1, Jan. 2004.
%
%Kindly report any suggestions or corrections to zhouwang@ieee.org
%
%----------------------------------------------------------------------
%
%Input : (1) img1: the first image being compared
%        (2) img2: the second image being compared
%        (3) roi: this is either a 4-vector with indexes
%            (xmin, xmax, ymin, ymax) defining a rectangular region of
%            interest or a binary mask (of the same size as the input
%            images) where the SSIM is computed. default is to compute
%            the SSIM on the full image
%        (4) K: constants in the SSIM index formula (see the above
%            reference). defualt value: K = [0.01 0.03]
%        (5) window: local window for statistics (see the above
%            reference). default widnow is Gaussian given by
%            window = fspecial('gaussian', 11, 1.5);
%        (6) L: dynamic range of the images. default: L = 255
%
%Output: (1) mssim: the mean SSIM index value between 2 images.
%            If one of the images being compared is regarded as 
%            perfect quality, then mssim can be considered as the
%            quality measure of the other image.
%            If img1 = img2, then mssim = 1.
%        (2) ssim_map: the SSIM index map of the test image. The map
%            has a smaller size than the input images. The actual size:
%            size(img1) - size(window) + 1.
%
%Default Usage:
%   Given 2 test images img1 and img2, whose dynamic range is 0-255
%
%   [mssim ssim_map] = ssim_index(img1, img2);
%
%Advanced Usage:
%   User defined parameters. For example
%
%   K = [0.05 0.05];
%   window = ones(8);
%   L = 100;
%   roi = [4, 7, 2, 10]; % or
%   roi = zeros(size(img1)); roi(4:7, 2:10) = 1; % with the same result
%   [mssim ssim_map] = ssim_index(img1, img2, roi, K, window, L);
%
%See the results:
%
%   mssim                        %Gives the mssim value
%   imshow(max(0, ssim_map).^4)  %Shows the SSIM index map
%
%========================================================================
%Modifications by Andreas Keil, Stanford University, 2011:
%- Enabled the selection of a region of interest where the similarity is
%  computed by introducing the new parameter roi.
%- Simplified argument checking.
%- Some other code simplifications/clarifications.



if nargin < 2
	error('At least two input images are needed!');
end

if ndims(img1) ~= 2 || ndims(img2) ~= 2
	error('Input images have to be 2-dimensional!');
end
[M N] = size(img1);
if ~all(size(img2) == [M N])
   error('Input images must have the same size!');
end

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

if nargin < 4
	K = [0.01 0.03];
end
if length(K) ~= 2 || K(1) < 0 || K(2) < 0
	error('K must be a 2-vector with non-negative numbers!');
end

if nargin < 5
	if ((M < 11) || (N < 11))
		error('Input images are too small for window!');
	end
	window = fspecial('gaussian', 11, 1.5);
end
[W1 W2] = size(window);
if (W1*W2) < 4 || (W1 > M) || (W2 > N)
	error('Window is either too small or larger than input images!');
end

if nargin < 6
	L = 255;                                  %
end

if nargin > 6
	error('Too many input parameters specified!');
end



C1 = (K(1)*L)^2;
C2 = (K(2)*L)^2;
window = window/sum(sum(window));
img1 = double(img1);
img2 = double(img2);

mu1   = filter2(window, img1, 'valid');
mu2   = filter2(window, img2, 'valid');
mu1_sq = mu1.^2;
mu2_sq = mu2.^2;
mu1_mu2 = mu1.*mu2;
sigma1_sq = filter2(window, img1.^2, 'valid') - mu1_sq;
sigma2_sq = filter2(window, img2.^2, 'valid') - mu2_sq;
sigma12 = filter2(window, img1.*img2, 'valid') - mu1_mu2;

numerator1 = 2*mu1_mu2 + C1;
numerator2 = 2*sigma12 + C2;
denominator1 = mu1_sq + mu2_sq + C1;
denominator2 = sigma1_sq + sigma2_sq + C2;

if C1 > 0 && C2 > 0
   ssim_map = (numerator1.*numerator2)./(denominator1.*denominator2);
else
   ssim_map = ones(size(mu1));
   index = (denominator1.*denominator2 > 0);
   ssim_map(index) = (numerator1(index).*numerator2(index))./(denominator1(index).*denominator2(index));
   index = (denominator1 ~= 0) & (denominator2 == 0);
   ssim_map(index) = numerator1(index)./denominator1(index);
end

% adapt ROI to filtered input images
crop1 = (W1-1)/2;
crop2 = (W2-1)/2;
[R1 R2] = size(roi);
roi = roi(1+crop1:R1-crop1, 1+crop2:R2-crop2);

% restrict computation of SSIM to ROI
mssim = mean(ssim_map(roi));
