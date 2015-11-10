function [sinoOut, geomOut] = upsamplingViewFeaturePreserving3( sinoIn, geomIn, usr, srchWidth, detDSR, threshold )
%% function [sinoOut, geomOut] = upsamplingViewFeaturePerserving( sinoIn,
% geomIn, usr, srchWidth, detDSR )
% Up sampling the sinogram in the azmuthal direction while perserving the
% small local structural features.
%
% Inputs:
%   sinoIn, geomIn, usr,
%   srchWidth - maximum search size in number of pixel
%   detDSR - detector down sampling rate
% Outputs:
%   sinoOut, geomOut
%
% Note: this is the 2D version for proof the concept
%
% Meng Wu at University of Erlangen Nuremberg
% 2014.10

[nv, nu, noViews ] = size(sinoIn);
geomOut = geomIn;

geomOut.noViews = (noViews-1)*usr;
sinoOut         = zeros( [nv, nu, geomOut.noViews], 'single');
geomOut.betas   = zeros(1, geomOut.noViews);
geomOut.couchZ  = zeros(1, geomOut.noViews);

featureShifts = zeros( [ nv, nu, noViews - 1], 'single');

fprintf('Featurs perserving interpolation: ')

for i = 0:noViews-2
    
    if mod(i, 20 ) == 0
        fprintf('(%i/%i)... ', i, noViews );
    end
    
        
    view1 = squeeze( sinoIn(:,:,i+1) );
    view2 = squeeze( sinoIn(:,:,i+2) );
    beta1 = geomIn.betas(i+1);
    beta2 = geomIn.betas(i+2);
    couchz1 = geomIn.couchZ(i+1);
    couchz2 = geomIn.couchZ(i+2);
    
    featureShifts( :, :, i + 1 ) = interpolateFeaturePerserving( view1, view2, srchWidth, detDSR, threshold );
    
    for j = 0:usr-1
        w1 = ( usr - j ) / usr ;
        w2 = j / usr;
        
        sinoOut( :, :, i*usr + j + 1 )    = w1 * view1 + w2 * view2 ;
        geomOut.betas( i*usr + j + 1  )   = w1 * beta1 + w2 * beta2 ;
        geomOut.couchZ( i*usr + j + 1  )  = w1 * couchz1 + w2 * couchz2;
    end
    
    
end


fprintf('done.\n')

sinoOut(:, :, end)    = sinoIn(:, :, end);
geomOut.betas(end)  = geomIn.betas(end);
geomOut.couchZ(end) = geomIn.couchZ(end);

% interpolation with featrue shifts
iu = (1:nu);
for i = 0 :  noViews-2
    
    view1 = sinoIn(:, :, i+1);
    view2 = sinoIn(:, :, i+2);
    
    for j = 1:usr-1
        
        view = sinoOut( :, :, i * usr + j + 1 );
        
        for iv = 1 : nv
            
            
            shift = featureShifts(iv, :, i+1);
            
            % interpoalte the two base slice
            slice1 = interp1( view1( iv, : ) , iu + shift );
            slice2 = interp1( view2( iv, : ) , iu - shift );
            
            c = shift * ( 2 * j / usr - 1 ); 
            
            % get the interpolated slice
            w1 = ( usr - j ) / usr ;
            w2 = j / usr;
            view( iv, : ) = interp1(  w1 *  slice1 +  w2 * slice2, iu + c );
            
        end
        
        sinoOut( :, :, i * usr + j + 1 ) = view;
        
    end
    
end


% for iv = 1 : size( sinoOut, 3 )
%     
%     center = sinoOut( end/2, end/2-8:end/2+7, iv );
%     
%     sinoOut(:,:,iv) = sinoOut(:,:,iv)  - mean( center(:) );
% end

end

%% function for feature perserving interpolation
function featureShifts = interpolateFeaturePerserving( view1, view2, srchWidth, detDSR, threshold )

compWidth = ceil( 12 / detDSR(1) );
srchWidth = ceil( srchWidth / detDSR(1) );

% spatial weights to drag sreaching maximum to the center
weights = exp( 0.5 * abs( - srchWidth : srchWidth ) /  srchWidth  );
featureShifts = zeros( size( view1 ) );

h = fspecial('Gaussian', [8 8], 2);
view1 = imfilter( view1, h , 'same' );
view2 = imfilter( view2, h , 'same' );

h2 = fspecial('Gaussian', [32 32], 8);
view1 = view1 - imfilter( view1, h2 , 'same' );
view2 = view2 - imfilter( view2, h2 , 'same' );

% down sample the pixel in the searching for timing consideration
subView1 = view1( ceil(detDSR(2)/2) : detDSR(2) : end , ceil(detDSR(1)/2) : detDSR(1) : end )';
subView2 = view2( ceil(detDSR(2)/2) : detDSR(2) : end , ceil(detDSR(1)/2) : detDSR(1) : end )';

subView = ( abs(subView1) + abs(subView2) ) / 2;
% subView1 = detectorPixelBinning( view1', [detDSR(1) detDSR(2)] );
% subView2 = detectorPixelBinning( view2', [detDSR(1) detDSR(2)] );

%subView1 = subView1 - mean( subView1(:) );
%subView2 = subView2 - mean( subView2(:) );


[nu, nv] = size( subView1 );

for j = 1 : nv 
    
    for i = round( nu / 8 )  + srchWidth + compWidth  :   round( nu * 7 / 8 )  - srchWidth - compWidth

        p0 = subView( i - compWidth : i + compWidth, j );
        
        scores = zeros( 1, 2*srchWidth + 1);
        for k = - srchWidth : srchWidth
            
            % find a segement in the opposite direction
            l1 = subView1( i - compWidth + k : i + compWidth + k, j );
            l2 = subView2( i - compWidth - k : i + compWidth - k, j );
                        
            % compute the match score
            scores( k + srchWidth +  1 ) = sum(  abs( l1(:) - l2(:) ) ) / sum( p0(:) ) ;
            
        end
        
        scores = scores .* weights;
        
        % find the weighted maximum
        [ ms, s ] = min( scores ) ;
        s = s - srchWidth - 1;
        
        if ms > threshold
            s = 0;
        end
       
        % record the interp featureShifts
        featureShifts( (j-1)*detDSR(2)+1 : j*detDSR(2) ,  (i-1)*detDSR(1)+1 : i*detDSR(1) ) = s * detDSR(1);
        
    end
    
end

% smooth a little bit
featureShifts = medfilt2( featureShifts, [3 5] );

h = fspecial('Gaussian', [8 8], 2);
featureShifts = imfilter( featureShifts, h , 'same' );


% h = figure(100);
% if ~ishandle(h)
% figure(100);
% end
% %quiver( featureShifts, zeros( size( featureShifts) ) );
% imdisp( featureShifts, [-srchWidth srchWidth] );
% pause(0.1 )


end

