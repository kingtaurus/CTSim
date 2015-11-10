function [ sinoOut, geomOut ] = truncationCorrectionEllipicalPatch( sinoIn, geomIn, mu, edgeWidth, rho, extendWidth, medfiltWidth )
% function [ sinoOut, geomOut ] = truncationCorrectionEllipicalPatch( sinoIn, geomIn, mu, edgeWidth, rho )
%
%   mu - attenuation coefficent of the base material ( usually water)
%   edgeWidth - number of pixel at the edge used for truncation detection
%   rho       - ratio of the axes ( default 1 )
%
%
% Meng Wu @ Stanford Univercity
% Created 2014.2
% Modified 2014.3 add 2d case

if nargin < 4
    edgeWidth = 3;
end

if nargin < 5
    rho = 1;
end

if nargin < 6
    extendWidth = round( geomIn.detSize(1) / 4 );
end

if nargin < 7
    medfiltWidth = 5;
end


a = rho * geomIn.detSpacing(1)  * mu / 10;

fittingWidth = edgeWidth * 2;
X = [ (1:fittingWidth)' ones(fittingWidth,1)];
Xinv = (X'*X)^(-1) * X';


if ndims( sinoIn ) == 3 %3D cases
    
    nu = geomIn.detSize(1);
    nv = geomIn.detSize(2);
    
    sinoOut = zeros( nv, nu + extendWidth * 2 , geomIn.noViews, 'single' );
    
    for iview = 1 : geomIn.noViews
        
        view = sinoIn(:,:,iview);
        %view = medfilt2(view,[medfiltWidth medfiltWidth], 'symmetric');
        
        h = fspecial('gaussian', 8, 2);
        view = imfilter(view, h, 'replicate' );
        
        viewExtrap = zeros( nv, nu + extendWidth * 2, 'single' );
        
        for iv = 1 : nv
            viewExtrap(iv,:) = extrapolatingSlice( view(iv,:), a, Xinv,  extendWidth, edgeWidth, fittingWidth  );
        end
        
        viewExtrap(:,extendWidth+1: extendWidth+nu ) = sinoIn(:,:,iview);
        
        if medfiltWidth > 1
            viewExtrap = medfilt2(viewExtrap,[medfiltWidth medfiltWidth]);
            % viewExtrap(:,extendWidth+edgeWidth+1: extendWidth+nu-edgeWidth ) = sinoIn(:,1+edgeWidth:end-edgeWidth,iview);
            viewExtrap(:,extendWidth+1: extendWidth+nu ) = sinoIn(:,:,iview);
        end
        
        sinoOut(:,:,iview) = viewExtrap;
        
    end
    
elseif ismatrix( sinoIn )  %2D cases
    
    nu = geomIn.detSize(1);
    h = fspecial('gaussian', 8, 2);
    sinoFilt = imfilter(sinoIn, h, 'replicate' );
    
    sinoOut = zeros( nu + extendWidth * 2 , geomIn.noViews, 'single' );
    
    for iview = 1 : geomIn.noViews
        sinoOut(:,iview) = extrapolatingSlice( sinoFilt(:,iview)', a, Xinv, extendWidth, edgeWidth, fittingWidth  );
        sinoOut(extendWidth+1: extendWidth+nu,: ) = sinoIn;
    end
    
    if medfiltWidth > 1
        sinoOut = medfilt2(sinoOut,[medfiltWidth 1]);
        sinoOut(extendWidth+edgeWidth+1: extendWidth+nu-edgeWidth,: ) = sinoIn(1+edgeWidth:end-edgeWidth,:);
    end
    
    
else
    fprintf('will be implemented in the future.\n')
    
end

geomOut = geomIn;
geomOut.detSize(1) = nu + extendWidth * 2;


end





function sliceExtrap = extrapolatingSlice( detRow, a, Xinv, extendWidth, edgeWidth, fittingWidth  )

nu = length( detRow );

sliceExtrap = zeros( 1, nu + extendWidth * 2, 'single' );

edgeValue = mean( detRow(1:edgeWidth) );

if  edgeValue > 0.05 % truncation happens at the left edge
    
    % number of extendent pixel if f
    t = round( edgeValue / a );
    if t < 3
        return;
    end
    
    % fitt the edge
    p = Xinv * detRow(1:fittingWidth)';
    edgeValue = p(2) + p(1);
    
    
    b = max( p(1), 0);
    if isnan(b) || isinf(b)
        b = 0;
    end
    b = min(b,  a * 2 );
    
    x = floor( t * a / (  a + b ) );
    x = min( x, extendWidth );
    
    ix = (-x:0);
    c = max(sqrt( t^2 - (t-x)^2 ), 1);
    sliceExtrap( extendWidth + ix  + 1 ) = edgeValue / c * sqrt( t^2 - (ix+x-t).^2 );
    
end


edgeValue =  mean( detRow(end-edgeWidth+1:end) );

if edgeValue > 0.05  % truncation happens at the right edge
    
    t = round( edgeValue / a );
    if t < 3
        return;
    end
    
    p = Xinv * detRow(end-fittingWidth+1:end)';
    edgeValue = p(2) + fittingWidth * p(1);
    
    b = max( - p(1), 0);
    if isnan(b) || isinf(b)
        b = 0;
    end
    b = min(b,  a * 2 );
    
    x = floor( t * a / ( a + b ) );
    x = min( x, extendWidth );
    ix = (0:x);
    c = max(sqrt( t^2 - (t-x)^2 ), 1);
    sliceExtrap( extendWidth+ nu + ix ) = edgeValue / c * sqrt( t^2 - (ix+t-x).^2 );
    
end

end

