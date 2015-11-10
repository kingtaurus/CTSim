function [mssim1, mssim2, rmse1, rmse2, stdTissue, R] = compare2GroundTruth( imgRecon, imgGt, roi, map, fileName )
% Compare reconstructio with ground truth. Two regions will be compared and
%   using SSIM and RMSE. STD of soft tissue will also be repoted
%  input:
%       imgRecon (in HU)
%       imgGt
%       roi
%       map
%   output:
%       mssim1, mssim2, rmse1, rmse2
%
% Meng Wu @ Stanford
% 2012


if nargin < 5
    fileName = 'Reconstruction';
end


if min(imgRecon(:)) > -800 || min(imgGt(:)) > -800
    fprintf('Warning: please use HU image. \n' );
end


if ndims( imgRecon ) == 3
    imgRecon = imgRecon( :,:, round( end/2 ) );
end


% other SSIM parameters
ssimK = [0.01 0.03];
ssimWindow = fspecial('gaussian', 11, 1.5);
ssimLKeV = 2000+1000;

% evaluate
mssim1 = my_ssim_index(imgRecon+1000, imgGt+1000, roi.ssimRoi1, ssimK, ssimWindow, ssimLKeV);
mssim2 = my_ssim_index(imgRecon+1000, imgGt+1000, roi.ssimRoi2, ssimK, ssimWindow, ssimLKeV);

rmse1 = computeRMSE(imgRecon, imgGt, roi.ssimRoi1);
rmse2 = computeRMSE(imgRecon, imgGt, roi.ssimRoi2);

R = 0;

d = abs(imgRecon(1:end-1, :) -  imgRecon(2:end, :));
R = R + sum( d( map.mapTissue(:)) );

d = abs(imgRecon(:, 1:end-1) -  imgRecon(:, 2:end));
R = R + sum( d( map.mapTissue(:)) );


fprintf('----------------------------------------------------------------------------\n');
fprintf('ROI 1 (betw. fillings): SSIM = %4.1f%%,  RMSE = %3i HU.\n', mssim1*100, round(rmse1));
fprintf('ROI 2 (outs. fillings): SSIM = %4.1f%%,  RMSE = %3i HU.\n', mssim2*100, round(rmse2));
if nargin > 3
    imgTissue = imgRecon( map.mapTissue(:) );
    stdTissue = std( imgTissue(:) );
    fprintf('Soft tissue           : STD  = %4.1fHU, TV = %4.3g.\n', stdTissue, R );
end
fprintf('----------------------------------------------------------------------------\n\n');

figure('Name', fileName , 'NumberTitle', 'off');
imshow( imgRecon', roi.windowHu);
% if nargin >= 5
%     saveas(gcf,  [ fileName '.eps'] );
% end

end