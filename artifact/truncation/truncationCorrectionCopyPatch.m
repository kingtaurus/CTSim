function [ sinoOut, geomOut ] = truncationCorrectionCopyPatch( sinoIn, geomIn, attSoft, edgeWidth, extendWidth )
% function [ sinoOut, geomOut ] = truncationCorrectionCopyPatch( sinoIn, geomIn, attenCoeff, edgeWidth, rho )
%
%   attenCoeff - attenuation coefficent of the base material ( usually water)
%   edgeWidth - number of pixel at the edge used for truncation detection
%   rho       - ratio of the axes ( in randon domain, round object is 2 )
%
%
% Meng Wu @ Stanford Univercity
% Created 2014.2
% Modified 2014.3 add 2d case


if nargin < 4
    edgeWidth = 3;
end


if nargin < 5
    extendWidth = round( geomIn.detSize(1) / 4 );
end


[ sinoOut, geomOut ] = truncationCorrectionWaterCylinderFitting( sinoIn, geomIn, attSoft, edgeWidth, extendWidth );

if ndims( sinoIn ) == 3
    
    nu = geomIn.detSize(1);
    nv = geomIn.detSize(2);
    
    for iview = 1 : geomIn.noViews
        
        view = sinoIn(:,:,iview);
        viewExtrap = sinoOut(:,:,iview);

        viewExtrap = medfilt2( viewExtrap, [3 5]);
        
        for iv = 1 : nv
            
            detRow = view(iv,:);
            detRowExtrap = viewExtrap(iv,:);
            
            edgeValueLeft = mean( detRow(1:edgeWidth) );
            edgeValueRight = mean( detRowExtrap(end-edgeWidth+1:end) );
            
            if  edgeValueLeft <  edgeValueRight &&  edgeValueRight > 0.05 % truncation happens at the right edge
                
                % find the most left pixel that has the same value as at edge
                % i =  find( detRow >= edgeValueRight, 1, 'first' );
                
                i = 1;
                peakValue = detRow(i);
                while i < nu && detRow(i+1) >= peakValue - 0.01 && peakValue  <= edgeValueRight || peakValue <= edgeValueRight * 0.5
                    i = i + 1;
                    peakValue = detRow(i);
                end
                
                scale = min( ( edgeValueRight / peakValue ), 2 );
                
                viewExtrap( iv, extendWidth + nu + 1 : end  ) = scale * flip( detRowExtrap( i + 1 :  i + extendWidth ) );
                
            elseif  edgeValueLeft > 0.05 % truncation happens at the left edge
                
                % find the most right pixel that has the same value as at edge
                % i =  find( detRow >= edgeValueLeft, 1, 'last' );
                
                i = nu;
                peakValue = detRow(i);
                while i > 1 && detRow(i-1) >= peakValue - 0.01 && peakValue <= edgeValueLeft  || peakValue <= edgeValueLeft * 0.5
                    i = i - 1;
                    peakValue = detRow(i);
                end
                
                scale = min( ( edgeValueLeft / peakValue ), 2 );
                
                viewExtrap( iv, 1:extendWidth ) = scale * flip(  viewExtrap( iv, i+extendWidth+1:i+2*extendWidth ) ) ;
                
            end
            
        end
        
        sinoOut(:,:,iview) = viewExtrap;
        sinoOut(:,extendWidth+1: extendWidth+nu, :) = sinoIn;
    end
    
    
elseif ismatrix( sinoIn ) %2D cases
    
    nu = geomIn.detSize(1);
    
    for iview = 1 : geomIn.noViews
        
        detRow = sinoIn(:,iview);
        detRowExtrap = sinoOut(:,iview);
        
        edgeValueLeft = mean( detRow(1:edgeWidth) );
        edgeValueRight = mean( detRowExtrap(end-edgeWidth+1:end) );
        
        if  edgeValueLeft <  edgeValueRight &&  edgeValueRight > 0.05 % truncation happens at the right edge
            
            % find the most left pixel that has the same value as at edge
            % i =  find( detRow >= edgeValueRight, 1, 'first' );
            i = 1;
            peakValue = detRow(i);
            while i < nu && detRow(i+1) >= peakValue - 0.01 && peakValue  <= edgeValueRight || peakValue <= edgeValueRight * 0.5
                i = i + 1;
                peakValue = detRow(i);
            end

            detRowExtrap( extendWidth + nu + 1 : end  ) = ( edgeValueRight / peakValue ) * flip( detRowExtrap( i + 1 :  i + extendWidth ) );
            
        elseif  edgeValueLeft > 0.05 % truncation happens at the left edge
            
            % find the most right pixel that has the same value as at edge
            % i =  find( detRow >= edgeValueLeft, 1, 'last' );
            i = nu;
            peakValue = detRow(i);
            while i > 1 && detRow(i-1) >= peakValue - 0.01 && peakValue <= edgeValueLeft  || peakValue <= edgeValueLeft * 0.5
                i = i - 1;
                peakValue = detRow(i);
            end

            detRowExtrap( 1:extendWidth ) = ( edgeValueRight / peakValue ) * flip(  detRowExtrap( i+extendWidth+1:i+2*extendWidth ) ) ;
            
        end
        
        sinoOut(:,iview) = detRowExtrap;
        
        
    end
    
else
    fprintf('will be implemented in the future.\n')
    
end

end


function A = flip( A )

A = A(:);
A = A(end:-1:1);
end
