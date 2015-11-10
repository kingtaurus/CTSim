function [ sinoOut, geomOut, truncatedView ] = truncationCorrectionWaterCylinderFitting( sinoIn, geomIn, mu, edgeWidth, extendWidth, truncatedView )
% function [ sinoOut, geomOut ] = truncationCorrectionWaterCylinderFitting( sinoIn, geomIn, mu, edgeWidth )
% Water cylinder fitting based truncation correction alogrithm for CT
% no total attenuation constrain version
%
% Based on:
% Hsieh, J., Chao, E., Thibault, J., Grekowicz, B., Horst, a., McOlash, S.,
% & Myers, T. J. (2004). A novel reconstruction algorithm to extend the CT
% scan field-of-view. Medical Physics, 31(9), 2385. doi:10.1118/1.1776673
%
% Meng Wu @ Stanford Univercity
% Created 2014

if nargin < 4
    edgeWidth = 3;
end

if nargin < 5
    extendWidth = round( geomIn.detSize(1) / 4 );
end

geomOut = geomIn;
geomOut.detSize(1) = geomIn.detSize(1) + extendWidth * 2;

if isfield( geomOut, 'PMat' )
    for iv = 1 : geomOut.noViews
        geomOut.PMat(:,:,iv) = [1 0 extendWidth; ...
            0  1 0; ...
            0  0  1 ] * squeeze( geomOut.PMat(:,:,iv) );
    end
end



% some parameters
threshold = 0.05;
mu = mu / 10 * geomIn.detSpacing(1) / 2;

% matrix used for fitting
fittingWidth = edgeWidth * 2;
X = [ (1:fittingWidth)' ones(fittingWidth,1)];
Xinv = (X'*X)^(-1) * X';

if nargin < 6
    truncatedView = detectTruncation( sinoIn, edgeWidth, threshold );
end

if ndims( sinoIn ) == 3 %3D cases
    
    % cosine weighting
    sinoIn = cosineRayWeighting( sinoIn, geomIn, 1 );
    
    nu = geomIn.detSize(1);
    nv = geomIn.detSize(2);
    
    sinoOut = zeros( nv, nu + extendWidth * 2 , geomIn.noViews, 'single' );
    
    for iview = 1 : geomIn.noViews
        
        if ~ truncatedView( iview )
            continue;
        end
        
        
        view = sinoIn(:,:,iview);
        viewExtrap = zeros( nv, nu + extendWidth * 2, 'single' );
        
        % [X, R] = cylinderFitting( mean( view ), mu, Xinv, edgeWidth, fittingWidth, threshold );
        
        hnv = round( nv / 4 );
        for iv = 1 : nv
            
            if mod( iv, 8 ) == 4 || iv == 1
                [X, R] = cylinderFitting( mean( view( max(1, iv-hnv): min(nv, iv+hnv), : )), mu, Xinv, edgeWidth, fittingWidth, threshold );
            end
            
            %viewExtrap(iv,:) = extrapolatingOneRow1( view(iv,:), mu, X,  extendWidth, edgeWidth, threshold  );
            viewExtrap(iv,:) = extrapolatingOneRow2( view(iv,:), X, R, extendWidth, edgeWidth, threshold  );
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
    sinoOut = cosineRayWeighting( sinoOut, geomOut, 0 );
    
elseif ismatrix( sinoIn )  %2D cases
    
    error('Not supported yet.');
    
end

end


function [X, R] = cylinderFitting( detRow, mu, Xinv, edgeWidth, fittingWidth, threshold )

X = [0 0];
R = [0 0];
if  mean( detRow(1:edgeWidth) ) > threshold % truncation happens at the left edge
    
    % fitt the edge
    p = Xinv * detRow(1:fittingWidth)';
    
    slope = p(1);
    P = p(1) + p(2);
    X(1) = ( slope * P ) / ( 4 * mu^2 );
    R(1) = sqrt(  P^2 / ( 4 * mu^2 ) + X(1)^2  );
end

if  mean( detRow(end-edgeWidth+1:end) ) > threshold  % truncation happens at the right edge
    
    % fitt the edge
    p = Xinv * detRow(end-fittingWidth+1:end)';
    
    slope = p(1);
    P = fittingWidth *  p(1) + p(2);
    X(2) = ( slope * P ) / ( 4 * mu^2 );
    R(2) = sqrt(  P^2 / ( 4 * mu^2 ) + X(2)^2  );
end

end


function sliceExtrap = extrapolatingOneRow1( detRow, mu, X, extendWidth, edgeWidth, threshold  )

nu = length( detRow );
sliceExtrap = zeros( 1, nu + extendWidth * 2, 'single' );

P = mean( detRow(1:edgeWidth) );
if  P > threshold % truncation happens at the left edge
    
    R = sqrt(  P^2 / ( 4 * mu^2 ) + X(1)^2  );
    x = abs( round(  R - X(1) ) );
    x = min( x, extendWidth );
    
    ix = ( - x + 1: 0 );
    d = abs( ix ) + X(1) ;
    
    sliceExtrap( extendWidth + ix ) = 2 * sqrt( R^2 - d.^2 ) * mu ;
    
end

P =  mean( detRow(end-edgeWidth+1:end) );
if P > threshold % truncation happens at the right edge
    
    R = sqrt(  P^2 / ( 4 * mu^2 ) + X(2)^2  );
    x = abs( round( ( R + X(2) )  ) );
    x = min( x, extendWidth );
    ix = ( 1 : x );
    d = ( ix  - 1 ) - X(2) ;
    
    sliceExtrap( extendWidth+ nu + ix ) = 2 * sqrt( R^2 - d.^2 ) * mu;
    
end

end

function sliceExtrap = extrapolatingOneRow2( detRow, X, R, extendWidth, edgeWidth, threshold  )

nu = length( detRow );
sliceExtrap = zeros( 1, nu + extendWidth * 2, 'single' );

P = mean( detRow(1:edgeWidth) );
if  abs(P) > threshold % truncation happens at the left edge
    
    x = abs( round(  R(1) - X(1) ) );
    x = min( x, extendWidth );
    
    ix = ( - x + 1: 0 );
    d = abs( ix ) + X(1) ;
    
    scale = P / sqrt( R(1)^2 - X(1)^2 );
    sliceExtrap( extendWidth + ix ) = scale * sqrt( R(1)^2 - d.^2 );
    
end

P =  mean( detRow(end-edgeWidth+1:end) );
if abs(P) > threshold % truncation happens at the right edge
    
    x = abs( round( ( R(2) + X(2) )  ) );
    x = min( x, extendWidth );
    ix = ( 1 : x );
    d = ( ix  - 1 ) - X(2) ;
    
    scale = P / sqrt( R(2)^2 - X(2)^2 );
    sliceExtrap( extendWidth+ nu + ix ) = scale * sqrt( R(2)^2 - d.^2 );
    
end

end



