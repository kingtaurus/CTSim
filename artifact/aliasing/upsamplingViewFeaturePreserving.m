function [sinoOut, geomOut] = upsamplingViewFeaturePreserving( sinoIn, geomIn, upSamplingRate, srchWidth, detDSR )
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


if round ( log2( upSamplingRate ) ) ~= log2( upSamplingRate )
    fprintf('Error: the up sampling rate must be power of 2. \n' );
end

% do twice up sampling recursively
while true
    [sinoOut, geomOut] = upsamplingViewBy2FeaturePerserving( sinoIn, geomIn, srchWidth, detDSR );
    
    upSamplingRate = upSamplingRate / 2;
    srchWidth = ceil( srchWidth / 2 );
    
    if upSamplingRate >= 2
        sinoIn = sinoOut;
        geomIn = geomOut;
       
    else
        break;
    end
    
end

end


%% function to do twice upsampling
function [sinoOut, geomOut] = upsamplingViewBy2FeaturePerserving( sinoIn, geomIn,  srchWidth, detDSR )

[nu, noViews ] = size(sinoIn);
geomOut = geomIn;

if geomIn.shortScan % short scan case
    
    geomOut.noViews = (noViews-1)*2+1;
    sinoOut         = zeros( [ nu, geomOut.noViews], 'single');
    geomOut.betas   = zeros(1, geomOut.noViews);
    geomOut.couchZ  = zeros(1, geomOut.noViews);
      
    shifts = zeros( [ nu, geomOut.noViews - noViews], 'single');
    
    for i = 0:noViews-2
        
        view1 = sinoIn(:,i+1);
        view2 = sinoIn(:,i+2);
        beta1 = geomIn.betas(i+1);
        beta2 = geomIn.betas(i+2);
        couchz1 = geomIn.couchZ(i+1);
        couchz2 = geomIn.couchZ(i+2);
        
        sinoOut( :, 2 * i + 1 ) = view1;
        sinoOut( :, 2 * i + 2 ) = ( view1 + view2 ) / 2;
        
        shifts( :, i + 1 ) = interpolateFeaturePerserving( view1, view2, srchWidth, detDSR );
        
        % just in the middle
        geomOut.betas( i*2 + 1 ) = beta1 ;
        geomOut.betas( i*2 + 2 ) = ( beta1 + beta2 ) / 2;
        geomOut.couchZ( i*2 + 1 ) = couchz1 ;
        geomOut.couchZ( i*2 + 2 ) = ( couchz1 + couchz2 ) / 2;
        
    end
    sinoOut(:,end)    = sinoIn(:,end);
    geomOut.betas(end)  = geomIn.betas(end);
    geomOut.couchZ(end) = geomIn.couchZ(end);
    
else % 360 scan case the first and last view are similar
    
    geomOut.noViews     = noViews*2;
    sinoOut             = zeros( [nu, geomOut.noViews], 'single');
    geomOut.betas       = zeros(1, geomOut.noViews);
    geomOut.couchZ      = zeros(1, geomOut.noViews);
    
    shifts = zeros( [ nu, geomOut.noViews - noViews], 'single');
    
    for i = 0 : noViews-1
        
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
        
        sinoOut( :, 2 * i + 1 ) = view1;
        sinoOut( :, 2 * i + 2 ) = ( view1 + view2 ) / 2;
        
       shifts( :, i + 1 ) = interpolateFeaturePerserving( view1, view2, srchWidth, detDSR );
       
        % just in the middle
        geomOut.betas( i*2 + 1 ) = beta1 ;
        geomOut.betas( i*2 + 2 ) = ( beta1 + beta2 ) / 2;
        geomOut.couchZ( i*2 + 1 ) = couchz1 ;
        geomOut.couchZ( i*2 + 2 ) = ( couchz1 + couchz2 ) / 2;
        
    end
    
    
end

shifts = medfilt2( shifts, [3 5]);


for i = 0 :  geomOut.noViews - noViews - 2 
    view1 = sinoIn(:,i+1);
    view2 = sinoIn(:,i+2);
    shift = shifts(:, i+1);
    
    sinoOut( :, 2 * i + 2 ) = ( view1( (1:nu)' + round( shift ) ) + view2( (1:nu)' - round( shift ) ) ) / 2;
end


end

%% function for feature perserving interpolation
function interpShifts = interpolateFeaturePerserving( view1, view2, srchWidth, detDSR, compSize, blurSize )

if nargin < 6
    compSize = 4;
    blurSize = 4;
end

srchWidth = ceil( srchWidth / detDSR  );


% linear interpolation, but will be changed later
view = ( view1 + view2 ) / 2;

% spatial weights to drag sreaching maximum to the center
weights = exp( - 2 * ( - srchWidth : srchWidth ).^2 /  srchWidth^2  );



% down sample the pixel in the searching for timing consideration
subView1 = detectorPixelBinning( view1, [detDSR 1] );
subView2 = detectorPixelBinning( view2, [detDSR 1] );
subView = ( subView1 + subView2 ) / 2;

L = length( view1 );
l = length( subView1 );

interpShifts = zeros( size( view1 ) );
for i = 1 + srchWidth + compSize  : l - srchWidth - compSize
    
    scores = zeros( 1, 2*srchWidth + 1);
    
    for k = - srchWidth : srchWidth
        
        % find a segement in the opposite direction
        l1 = subView1( i - compSize + k : i + compSize + k  );
        l2 = subView2( i - compSize - k : i + compSize - k  );
        
        %l1 = l1 - mean( l1 );
        %l2 = l2 - mean( l2 );
        
        % compute the match score
        scores( k + srchWidth +  1 ) = sum( l1 .* l2) ;
        
    end
    
    % find the weighted maximum
    [~,j] = max( scores .* weights );
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
h = fspecial('Gaussian', [blurSize 1], blurSize/2 );
interpShifts = conv( interpShifts, h , 'same' );



end

