function [sinoOut, geomOut] = upsamplingViews( sinoIn, geomIn, ratio )
% [sinoOut, geomOut] = upsamplingViews( sinoIn, geomIn, ratio )
% add interpolated views using linear interpolations
%
% Meng Wu @ Stanford University
% Created 2013
% Modified 2014.3 to add couch position z



if ndims( sinoIn ) == 3
    [nv, nu, noViews ] = size(sinoIn);
    geomOut = geomIn;
    if geomIn.shortScan % short scan case
        
        geomOut.noViews = (noViews-1)*ratio+1;
        sinoOut         = zeros( [nv, nu, geomOut.noViews], 'single');
        geomOut.betas   = zeros(1, geomOut.noViews);
        geomOut.couchZ  = zeros(1, geomOut.noViews);
        
        for i = 0:noViews-2
            
            view1 = sinoIn(:,:,i+1);
            view2 = sinoIn(:,:,i+2);
            beta1 = geomIn.betas(i+1);
            beta2 = geomIn.betas(i+2);
            couchz1 = geomIn.couchZ(i+1);
            couchz2 = geomIn.couchZ(i+2);
            
            for j = 0:ratio-1
                sinoOut( :,:, i*ratio + j + 1 )     = ( ratio - j ) / ratio * view1 + j / ratio * view2 ;
                geomOut.betas( i*ratio + j + 1  )   = ( ratio - j ) / ratio * beta1 + j / ratio * beta2 ;
                geomOut.couchZ( i*ratio + j + 1  )  = ( ratio - j ) / ratio * couchz1 + j / ratio * couchz2;
            end
            
        end
        sinoOut(:,:,end)    = sinoIn(:,:,end);
        geomOut.betas(end)  = geomIn.betas(end);
        geomOut.couchZ(end) = geomIn.couchZ(end);
        
    else % 360 scan case the first and last view are similar
        
        geomOut.noViews     = noViews*ratio;
        sinoOut             = zeros( [nv, nu, geomOut.noViews], 'single');
        geomOut.betas       = zeros(1, geomOut.noViews);
        geomOut.couchZ      = zeros(1, geomOut.noViews);
        
        for i = 0:noViews-1
            
            view1 = sinoIn(:,:,i+1);
            beta1 = geomIn.betas(i+1);
            couchz1 = geomIn.couchZ(i+1);
            
            if i == noViews-1
                beta2 = geomIn.betas(1);
                view2 = sinoIn(:,:,1);
                couchz2 = geomIn.couchZ(1);
                if beta2 > beta1
                    beta1 = beta1 + 2 * pi;
                else
                    beta2 = beta2 + 2 * pi;
                end
            else
                view2 = sinoIn(:,:,i+2);
                beta2 = geomIn.betas(i+2);
                couchz2 = geomIn.couchZ(i+2);
            end
            
            for j = 0:ratio-1
                sinoOut( :,:, i*ratio + j + 1 )     = ( ratio - j ) / ratio * view1 + j / ratio * view2 ;
                geomOut.betas( i*ratio + j + 1  )   = ( ratio - j ) / ratio * beta1 + j / ratio * beta2 ;
                geomOut.couchZ( i*ratio + j + 1  )  = ( ratio - j ) / ratio * couchz1 + j / ratio * couchz2;
            end
            
        end
        
    end
    
else
    
    [nu, noViews ] = size(sinoIn);
    geomOut = geomIn;
    
    if geomIn.shortScan % short scan case
        
        geomOut.noViews = (noViews-1)*ratio+1;
        sinoOut         = zeros( [ nu, geomOut.noViews], 'single');
        geomOut.betas   = zeros(1, geomOut.noViews);
        geomOut.couchZ  = zeros(1, geomOut.noViews);
        
        for i = 0:noViews-2
            
            view1 = sinoIn(:,i+1);
            view2 = sinoIn(:,i+2);
            beta1 = geomIn.betas(i+1);
            beta2 = geomIn.betas(i+2);
            couchz1 = geomIn.couchZ(i+1);
            couchz2 = geomIn.couchZ(i+2);
            
            for j = 0:ratio-1
                sinoOut( :, i*ratio + j + 1 )     = ( ratio - j ) / ratio * view1 + j / ratio * view2 ;
                geomOut.betas( i*ratio + j + 1  )   = ( ratio - j ) / ratio * beta1 + j / ratio * beta2 ;
                geomOut.couchZ( i*ratio + j + 1  )  = ( ratio - j ) / ratio * couchz1 + j / ratio * couchz2;
            end
            
        end
        sinoOut(:,end)    = sinoIn(:,end);
        geomOut.betas(end)  = geomIn.betas(end);
        geomOut.couchZ(end) = geomIn.couchZ(end);
        
    else % 360 scan case the first and last view are similar
        
        geomOut.noViews     = noViews*ratio;
        sinoOut             = zeros( [nu, geomOut.noViews], 'single');
        geomOut.betas       = zeros(1, geomOut.noViews);
        geomOut.couchZ      = zeros(1, geomOut.noViews);
        
        for i = 0:noViews-1
            
            view1 = sinoIn(:,i+1);
            beta1 = geomIn.betas(i+1);
            couchz1 = geomIn.couchZ(i+1);
            
            if i == noViews-1
                beta2 = geomIn.betas(1);
                view2 = sinoIn(:,1);
                couchz2 = geomIn.couchZ(1);
                if beta2 > beta1
                    beta1 = beta1 + 2 * pi;
                else
                    beta2 = beta2 + 2 * pi;
                end
            else
                view2 = sinoIn(:,i+2);
                beta2 = geomIn.betas(i+2);
                couchz2 = geomIn.couchZ(i+2);
            end
            
            for j = 0:ratio-1
                sinoOut( :, i*ratio + j + 1 )     = ( ratio - j ) / ratio * view1 + j / ratio * view2 ;
                geomOut.betas( i*ratio + j + 1  )   = ( ratio - j ) / ratio * beta1 + j / ratio * beta2 ;
                geomOut.couchZ( i*ratio + j + 1  )  = ( ratio - j ) / ratio * couchz1 + j / ratio * couchz2;
            end
            
        end
        
    end
    
    
    
end



