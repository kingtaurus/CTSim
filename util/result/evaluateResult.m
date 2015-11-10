function [rmseAll, tvAll, ssimAll, rmseTissue, tvTissue, stdTissue] = evaluateResult( imgAtt, imgGt, spectrum, map, methodName, outputDir, show )
% Evaluted image reconstruction results
%
% Meng Wu
% 2012 - 2013

if nargin < 7
    show = 0;
end


if strcmp( methodName, 'Gt')
    imgHu = imgAtt;
    %roi.windowHu = [0 7];
else
    
    if ndims( imgAtt ) == 3
        imgCenter = squeeze( imgAtt(:,:, ceil(size(imgAtt,3)/2)+1 ) );
        
        imgHuTop = convertMonoAttToHu( squeeze( imgAtt(:,:,1 ) ), spectrum);
        imgHuBottom = convertMonoAttToHu( squeeze( imgAtt(:,:,end ) ), spectrum);
    else
        imgCenter = imgAtt;
    end
    
    imgCenter =  imgMonoAttCorrection( imgCenter, spectrum, map );
    imgHu = convertMonoAttToHu( imgCenter, spectrum);
    
    % root mean square error
    imgTissue = imgHu( map.mapTissue(:) );
    stdTissue = std( imgTissue(:) );
    rmseAll  = std( imgHu(:) - imgGt(:) );
    rmseTissue = std( imgHu( map.mapTissue(:) ) - imgGt( map.mapTissue(:) ) );
    
    % other SSIM parameters
    ssimK = [0.01 0.03];
    ssimWindow = fspecial('gaussian', 11, 1.5);
    ssimLKeV = 2000+1000;
    ssimRoi = [1 size(imgHu, 1) 1 size(imgHu, 2)];
    ssimAll = my_ssim_index(imgHu+1000, imgGt+1000, ssimRoi, ssimK, ssimWindow, ssimLKeV);
    
    % total variation
    tvTissue   = 0;
    tvAll      = 0;
    
    d = abs(imgHu(1:end-1, :) -  imgHu(2:end, :));
    mapTissue = map.mapTissue(1:end-1, :);
    tvTissue = tvTissue + sum( d( mapTissue(:)) );
    tvAll = tvAll + sum( d(:) );
    
    d = abs(imgHu(:, 1:end-1) -  imgHu(:, 2:end));
    mapTissue = map.mapTissue(:,1:end-1);
    tvTissue = tvTissue + sum( d( mapTissue(:)) );
    tvAll = tvAll + sum( d(:) );
    
    fprintf('----------------------------------------------------------------------------\n');
    fprintf('Enntire img           : RMSE = %3.2fHU, \tTV = %4.3g, \tSSIM = %4.1f.\n', rmseAll, tvAll, ssimAll * 100 );
    fprintf('Soft tissue           : RMSE = %3.2fHU, \tTV = %4.3g, \tSTD  = %4.1fHU.\n', rmseTissue, tvTissue, stdTissue );
    fprintf('----------------------------------------------------------------------------\n\n');
    
    
end

if show
    
    figure('Name', methodName, 'NumberTitle', 'off');
    imshow( imgHu', [-400 800]);
    saveas(gcf, [outputDir methodName '.eps']);
    
    if show == 2
        figure('Name', [ methodName '_Top'], 'NumberTitle', 'off');
        imshow( imgHuTop', [-400 800]);
        saveas(gcf, [outputDir methodName '_top.eps']);
        
        figure('Name', [ methodName '_Bottom'], 'NumberTitle', 'off');
        imshow( imgHuBottom', [-400 800]);
        saveas(gcf, [outputDir methodName '_bottom.eps']);
    end

end

pause(2);
