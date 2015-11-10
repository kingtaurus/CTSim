function [ sinoOut, geomOut, truncatedView ] = truncationCorrectionReferenceDetector( sinoIn, geomIn, sinoRefLeft, sinoRefRight, edgeWidth, extendWidth )
% function [ sinoOut, geomOut ] = truncationCorrectionExtendedDetector( sinoIn, geomIn, mu, edgeWidth )
% Water cylinder fitting based truncation correction alogrithm for CT
% no total attenuation constrain version
%
%
% Meng Wu @ Stanford Univercity
% Created 2014


geomOut = geomIn;
geomOut.detSize(1) = geomIn.detSize(1) + extendWidth * 2;

% some parameters
threshold = 0.05;

if ndims( sinoIn ) == 3 %3D cases
    
    % cosine weighting
    % sinoIn = cosineRayWeighting( sinoIn, geomIn, 1 );
    
    nu = geomIn.detSize(1);
    nv = geomIn.detSize(2);
    
    sinoOut = zeros( nv, nu + extendWidth * 2 , geomIn.noViews, 'single' );
    
    for iview = 1 : geomIn.noViews
        
        view = sinoIn(:,:,iview);
        viewExtrap = zeros( nv, nu + extendWidth * 2, 'single' );
        
        for iv = 1 : nv
            viewExtrap(iv,:) = interpolateOneRow( view(iv,:), ...
                sinoRefLeft(:,iview), sinoRefRight(:,iview), extendWidth, edgeWidth, threshold );
        end
        viewExtrap(:,extendWidth+1: extendWidth+nu) = view;
        
        % smooth the edge
        viewExtrap(:,extendWidth+1: extendWidth+nu ) = view;
        h = fspecial('gaussian', 8, 4);
        viewExtrap = imfilter(viewExtrap, h, 'replicate' );
        
        sinoOut(:,:,iview) = viewExtrap;
        % cosine weighting
    end
    
    sinoOut(:,extendWidth+1: extendWidth+nu, :) = sinoIn;
    % inverse cosine weighting
    % sinoOut = cosineRayWeighting( sinoOut, geomOut, 0 );
    
elseif ismatrix( sinoIn )  %2D cases
    
    error('Not supported yet.');
    
end

end


function sliceExtrap = interpolateOneRow( detRow, detRowLeft, detRowRight, extendWidth, edgeWidth, threshold  )

nu = length( detRow );
sliceExtrap = zeros( 1, nu + extendWidth * 2, 'single' );

P = mean( detRow(1:edgeWidth) );
E = mean( detRowLeft( end - edgeWidth + 1 : end ) );

if E > threshold && P > threshold
    scale = P / E ;
else
    scale = 1;
end

scale = min( scale, 1.5 );
scale = max( scale, 0.8 );

sliceExtrap( 1:extendWidth ) = scale * detRowLeft  ;


P = mean( detRow(end-edgeWidth+1:end) );
E = mean( detRowRight( 1 : edgeWidth ));

if E > threshold && P > threshold
    scale = P / E ;
else
    scale = 1;
end

scale = min( scale, 1.5 );
scale = max( scale, 0.8 );

sliceExtrap( end - extendWidth + 1 : end ) = scale * detRowRight;


end



