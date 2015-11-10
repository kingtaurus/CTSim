function [truncatedView, totalAttenuations] = detectTruncation( sinoIn, edgeWidth, threshold )
% function [truncatedView, totalAttenuations] = detectTruncation( sinoIn, edgeWidth, threshold )
% Detecting truncation projection in sinogram
%
% Meng Wu at Stanford University
% 2014

if ndims( sinoIn ) == 3 %3D cases
    
    noViews = size( sinoIn, 3 );
    totalAttenuations = zeros( [noViews, 1] );
    truncatedView = zeros( [noViews, 1] );
    
    for iview = 1 : noViews
        
        
        view = squeeze( sinoIn(:,:,iview) );
        totalAttenuations( iview ) = mean( sum( view, 2 ) );
        
        % looks at two edege
        maxEdgeLeft     = max( mean( view(:, 1:edgeWidth ),  2) );
        maxEdgeRight    = max( mean( view(:, end-edgeWidth+1 : end ),  2) );
        
        % recode the trucation for that projection
        if maxEdgeLeft > threshold || maxEdgeRight > threshold
            truncatedView( iview ) = 1;
        end
        
        
    end
    
elseif ismatrix( sinoIn )  %2D cases

    noViews = size( sinoIn, 2 );
    totalAttenuations = zeros( [noViews, 1] );
    truncatedView = zeros( [noViews, 1] );
    
    for iview = 1 : noViews
        
        view = squeeze( sinoIn(:,iview) );
        totalAttenuations( iview ) = median( sum( view ) );
        
        % looks at two edege
        maxEdgeLeft     = mean( view(1:edgeWidth ) );
        maxEdgeRight    = mean( view(end-edgeWidth+1 : end ) );
        
        % recode the trucation for that projection
        if maxEdgeLeft > threshold || maxEdgeRight > threshold
            truncatedView( iview ) = 1;
        end        
    end

end

end