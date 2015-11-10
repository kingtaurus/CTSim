function [imgHu, mssim1, mssim2, rmse1, rmse2, softSTD, softTV] = evaluateResultKvmv( imgAtt, imgGt, spectrum, map, roi, methodName, outputDir, show )
% Evaluted image reconstruction results
%
% Meng Wu
% 2012 - 2013

if nargin < 8
    show = 0;
end


if strcmp( methodName, 'Gt')
    imgHu = imgAtt;
    %roi.windowHu = [0 7];
else
    
    if length( size(imgAtt) ) == 3
        imgAtt = squeeze( imgAtt(:,:, ceil(size(imgAtt,3)/2) ) );
    end
    
    imgAtt =  imgMonoAttCorrection( imgAtt, spectrum, map );
    imgHu = convertMonoAttToHu( imgAtt, spectrum);
    
    [mssim1, mssim2, rmse1, rmse2, softSTD, softTV] = compare2GroundTruth( imgHu, imgGt, roi, map );
end

if show
    
    figure('Name', methodName, 'NumberTitle', 'off');
    imshow( imgHu', roi.windowHu);
    saveas(gcf, [outputDir methodName '.eps']);
    
    figure('Name', [methodName ' ROI1'], 'NumberTitle', 'off');
    imgROI = imgHu(roi.ssimRoi1(1): roi.ssimRoi1(2), roi.ssimRoi1(3):roi.ssimRoi1(4));
    imgROI = imgROI - mean( imgROI(:));
    imshow( imgROI' , [-200 200]);
    saveas(gcf, [outputDir methodName '_ROI1.eps']);

    figure('Name', [methodName ' ROI2'], 'NumberTitle', 'off');
    imgROI = imgHu(roi.ssimRoi2(1): roi.ssimRoi2(2), roi.ssimRoi2(3):roi.ssimRoi2(4));
    imgROI = imgROI - mean( imgROI(:));
    imshow( imgROI' , [-200 200]);
    saveas(gcf, [outputDir methodName '_ROI2.eps']);
end

pause(2);
