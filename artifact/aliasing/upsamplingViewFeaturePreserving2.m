function [sinoOut, geomOut] = upsamplingViewFeaturePreserving2( sinoIn, geomIn, upSamplingRate, srchWidth, detDSR )
%% function [sinoOut, geomOut] = upsamplingViewFeaturePerserving( sinoIn,
% geomIn, upSamplingRate, srchWidth, detDSR )
% Up sampling the sinogram in the azmuthal direction while perserving the
% small local structural features.
%
% Inputs:
%   sinoIn, geomIn, upSamplingRate,
%   srchWidth - maximum search size in number of pixel
%   detDSR - detector down sampling rate
% Outputs:
%   sinoOut, geomOut
%
% Note: this is the 2D version for proof the concept
%
% Meng Wu at University of Erlangen Nuremberg
% 2014.10

[nu, noViews ] = size(sinoIn);
geomOut = geomIn;


geomOut.noViews = (noViews-1)*upSamplingRate;
sinoOut         = zeros( [ nu, geomOut.noViews], 'single');
geomOut.betas   = zeros(1, geomOut.noViews);
geomOut.couchZ  = zeros(1, geomOut.noViews);

shifts = zeros( [ nu,  noViews - 1], 'single');

for i = 0:noViews-2
    
    view1 = sinoIn(:,i+1);
    view2 = sinoIn(:,i+2);
    beta1 = geomIn.betas(i+1);
    beta2 = geomIn.betas(i+2);
    couchz1 = geomIn.couchZ(i+1);
    couchz2 = geomIn.couchZ(i+2);
    
    shifts( :, i + 1 ) = interpolateFeaturePerserving( view1, view2, srchWidth, detDSR );
    
    
    for j = 0:upSamplingRate-1
        sinoOut( :, i*upSamplingRate + j + 1 )     = ( upSamplingRate - j ) / upSamplingRate * view1 + j / upSamplingRate * view2 ;
        geomOut.betas( i*upSamplingRate + j + 1  )   = ( upSamplingRate - j ) / upSamplingRate * beta1 + j / upSamplingRate * beta2 ;
        geomOut.couchZ( i*upSamplingRate + j + 1  )  = ( upSamplingRate - j ) / upSamplingRate * couchz1 + j / upSamplingRate * couchz2;
    end
    
    
end
sinoOut(:,end)    = sinoIn(:,end);
geomOut.betas(end)  = geomIn.betas(end);
geomOut.couchZ(end) = geomIn.couchZ(end);

shifts = medfilt2( shifts, [upSamplingRate upSamplingRate]);
h = fspecial('Gaussian', [upSamplingRate upSamplingRate], upSamplingRate/2 );
shifts = imfilter( shifts, h , 'same' );



for i = 0 :  noViews-2
    view1 = sinoIn(:,i+1);
    view2 = sinoIn(:,i+2);
    shift = shifts(:, i+1);
    
    for j = 1:upSamplingRate-1
        sinoOut( (1:nu)' - round( shift * ( 2 * j /upSamplingRate - 1) ), i * upSamplingRate + j + 1 ) = ...
            ( upSamplingRate - j ) / upSamplingRate * view1( (1:nu)' + round( shift ) ) ...
            + j / upSamplingRate * view2( (1:nu)' -  round( shift ) ) ;
        
    end
    
 end

%sinoOut = medfilt2(sinoOut, [3 3]);
end

%% function for feature perserving interpolation
function interpShifts = interpolateFeaturePerserving( view1, view2, srchWidth, detDSR, compSize, blurSize )

if nargin < 6
    compSize = 6;
    blurSize = 4;
end

srchWidth = ceil( srchWidth / detDSR  );

% spatial weights to drag sreaching maximum to the center
weights = exp( - 2 * ( - srchWidth : srchWidth ).^2 /  srchWidth^2  );

% down sample the pixel in the searching for timing consideration
subView1 = detectorPixelBinning( view1, [detDSR 1] );
subView2 = detectorPixelBinning( view2, [detDSR 1] );

% linear interpolation, but will be changed later
subView = ( subView1 + subView2 ) / 2;
l = length( subView1 );

interpShifts = zeros( size( view1 ) );
for i = 1 + srchWidth + compSize  :  l - srchWidth - compSize 
       
    if i > round( 0.4 * l ) && i < round( 0.6 * l )
        continue;
    end
    
    scores = zeros( 1, 2*srchWidth + 1);
    for k = - srchWidth : srchWidth
        
        % find a segement in the opposite direction
        l1 = subView1( i - compSize + k : i + compSize + k  );
        l2 = subView2( i - compSize - k : i + compSize - k  );
        
%         l1 = l1 - mean( l1 );
%         l2 = l2 - mean( l2 );

        % compute the match score
        scores( k + srchWidth +  1 ) = sum( l1 .* l2);
        
    end

    scores = scores .* weights;
    
    % find the weighted maximum
    [ ~ ,j] = max( scores );
    j = j - srchWidth - 1;
    
    % record the interp shifts
    interpShifts( (i-1) * detDSR + 1: i * detDSR ) = j * detDSR;
    
    % interpolate the subsampled view for fun
    subView( i ) = ( subView1(i + j ) + subView2(i - j ) ) / 2;
    
end

% % blur the probjection a little bit for denoising
% h = fspecial('Gaussian', [blurSize 1], blurSize/2 );
% view1 = conv( view1, h, 'same' ) ;
% view2 = conv( view2, h, 'same' ) ;


% smooth a little bit



end

