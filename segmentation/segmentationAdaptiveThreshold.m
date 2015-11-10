function [segThr, segCtr, segStd] = segmentationAdaptiveThreshold( img, segThr)


nitn = 4;
numClusters     = length( segThr ) + 1;
segCtr    = zeros( 1, numClusters );
segStd     = zeros( 1, numClusters );

for itn = 1:nitn
    
    for ic = 1:numClusters
        
        % thresholding
        if ic == 1
            map = img < segThr(ic);
        elseif ic == numClusters
            map = img > segThr(ic-1);
        else
            map = img > segThr(ic-1) & img < segThr(ic);
        end
        
        % find cluster center and standard deviation
        segCtr(ic)    = mean(  img( map(:) ) );
        segStd(ic)     = std(  img( map(:) ) );
    end
    
    % update the thresholds
    for ic = 1:numClusters-1
        
        a = segStd(ic);
        b = segStd(ic+1);
        
        alpha = segCtr(ic);
        beta = segCtr(ic+1);
        
        % find a threshold
        for k = 3 : -1 : 0
            
            aa = alpha + k * a ;
            bb = beta - k * b;
            
            if aa < bb
                segThr(ic) = ( aa + bb ) / 2;
                break;
            end
        end
        
    end
    
end
segThr( isnan(segThr)) = max( segThr );
segCtr( isnan(segCtr)) = max( segCtr );

end