function [ sinoOut, geomOut, truncatedView ] = truncationCorrectionWaterCylinderFittingTotalConsistency( sinoIn, geomIn, mu, edgeWidth, extendWidth )
% function [ sinoOut, geomOut ] = truncationCorrectionWaterCylinderFittingTotalConsistency( sinoIn, geomIn, mu, edgeWidth )
% Water cylinder fitting based truncation correction alogrithm for CT
% using total attenuation constrain scaling version
%
% Based on:
% Hsieh, J., Chao, E., Thibault, J., Grekowicz, B., Horst, a., McOlash, S.,
% & Myers, T. J. (2004). A novel reconstruction algorithm to extend the CT
% scan field-of-view. Medical Physics, 31(9), 2385. doi:10.1118/1.1776673
% Meng Wu @ Stanford Univercity
%
% Modified 2014

if nargin < 4
    edgeWidth = 3;
end

if nargin < 5
    extendWidth = round( geomIn.detSize(1) / 4 );
end

geomOut = geomIn;
geomOut.detSize(1) = geomIn.detSize(1) + extendWidth * 2;

% some parameters
threshold = 0.05;
mu = mu / 10 * geomIn.detSpacing(1) / 2;

% matrix used for fitting
fittingWidth = edgeWidth * 2;
X = [ (1:fittingWidth)' ones(fittingWidth,1)];
Xinv = (X'*X)^(-1) * X';

if ndims( sinoIn ) == 3 %3D cases
    
    % cosine weighting
    sinoIn = cosineRayWeighting( sinoIn, geomIn, 1 );
    
    nu = geomIn.detSize(1);
    nv = geomIn.detSize(2);
    
    [truncatedView, totalAtt] = detectTruncation( sinoIn, edgeWidth, threshold );
    totalAttDiff = fillTotalAttenuation( truncatedView,  totalAtt ) - totalAtt;
    
    %figure; plot( totalAtt ); hold on; plot( totalAtt + totalAttDiff, ':' ); axis([1 geomIn.noViews 500 1000]);
    
    sinoOut = zeros( nv, nu + extendWidth * 2 , geomIn.noViews, 'single' );
    
    for iview = 1 : geomIn.noViews
        
        if ~ truncatedView( iview )
            continue;
        end
        
        view = sinoIn(:,:,iview);
        viewExtrap = zeros( nv, nu + extendWidth * 2, 'single' );
        [X, R, Eta] = cylinderFitting( mean( view ), mu, Xinv, edgeWidth, fittingWidth, threshold, totalAttDiff(iview) );
        
        hnv = round( nv / 4 );
        for iv = 1 : nv
            
            if mod( iv, 8 ) == 4 || iv == 1
                [X, R, Eta] = cylinderFitting( mean( view( max(1, iv-hnv): min(nv, iv+hnv), : )), mu, Xinv, edgeWidth, fittingWidth, threshold, totalAttDiff(iview) );
            end
            viewExtrap(iv,:) = extrapolatingOneRow( view(iv,:), X, R, Eta, extendWidth, edgeWidth, threshold  );
        end
        viewExtrap(:,extendWidth+1: extendWidth+nu) = view;
        
        % smooth the edge
        viewExtrap(:,extendWidth+1: extendWidth+nu ) = view;
        h = fspecial('gaussian', 8, 4);
        viewExtrap = imfilter(viewExtrap, h, 'replicate' );
        
        sinoOut(:,:,iview) = viewExtrap;
        
    end
    
    sinoOut(:,extendWidth+1: extendWidth+nu, :) = sinoIn;
    % inverse cosine weighting
    sinoOut = cosineRayWeighting( sinoOut, geomOut, 0 );
    
    
    
elseif ismatrix( sinoIn )  %2D cases
    
    disp('Not supported yet.');
    
end

end



function [X, R, Eta] = cylinderFitting( detRow, mu, Xinv, edgeWidth, fittingWidth, threshold, attDiff )

X = [0 0];
R = [0 0];

totalPatch = 0;

if  mean( detRow(1:edgeWidth) ) > threshold % truncation happens at the left edge
    
    % fitt the edge
    p = Xinv * detRow(1:fittingWidth)';
    
    slope = p(1);
    P = p(1) + p(2);
    X(1) = ( slope * P ) / ( 4 * mu^2 );
    R(1) = sqrt(  P^2 / ( 4 * mu^2 ) + X(1)^2  );
    
    d = abs( - abs( round(  R(1) - X(1) ) ) + 1 : 0 ) + X(1) ;
    totalPatch = totalPatch + sum( 2 * mu * sqrt( R(1)^2 - d.^2 ) );
    
end

if  mean( detRow(end-edgeWidth+1:end) ) > threshold  % truncation happens at the right edge
    
    % fitt the edge
    p = Xinv * detRow(end-fittingWidth+1:end)';
    
    slope = p(1);
    P = fittingWidth *  p(1) + p(2);
    X(2) = ( slope * P ) / ( 4 * mu^2 );
    R(2) = sqrt(  P^2 / ( 4 * mu^2 ) + X(2)^2  );
    
    d = ( 0 : abs( round( ( R(2) + X(2) ) ) ) - 1 ) - X(2) ;
    totalPatch = totalPatch + sum( 2 * mu * sqrt( R(2)^2 - d.^2 ) );
    
end

Eta = attDiff / totalPatch;

Eta = min( Eta, 5 );
Eta = max( Eta, 0.2 );

end


function sliceExtrap = extrapolatingOneRow( detRow, X, R, Eta, extendWidth, edgeWidth, threshold  )

nu = length( detRow );
sliceExtrap = zeros( 1, nu + extendWidth * 2, 'single' );

P = mean( detRow(1:edgeWidth) );
if  P > threshold % truncation happens at the left edge
    
    x = round( abs(   R(1) - X(1) ) * Eta )  ;
    x = min( x , extendWidth );
    
    ix = ( - x + 1: 0 );
    d = abs( ix ) / Eta + X(1) ;
    
    scale = P / sqrt( R(1)^2 - X(1)^2 );
    sliceExtrap( extendWidth + ix ) = scale * sqrt( R(1)^2 - d.^2 );
    
end

P =  mean( detRow(end-edgeWidth+1:end) );
if P > threshold % truncation happens at the right edge
    
    x = round( abs( ( R(2) + X(2) ) ) * Eta );
    x = min( x , extendWidth );
    ix = ( 1 : x );
    d = ( ix  - 1 )  / Eta - X(2) ;
    
    scale = P / sqrt( R(2)^2 - X(2)^2 );
    sliceExtrap( extendWidth+ nu + ix ) = scale * sqrt( R(2)^2 - d.^2 );
    
end

end



