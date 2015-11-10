function [xs, phis, regs] = reconPwlsApproxPathSeekingOld( y, w, geom, pfun, delta, numos, img0, img1, dv, updateRate, noFrames )
% Penalized-Weight-Least-Squares recosntruction for single energy sinograms
% using Nonuniform SQS algorithm
% input:
%       y       - log sinogram
%       w       - weights for pwls
%       geom    - system geometry
%       pfun    - roughness penalty function (default 'huber')
%       delta   - parameter for penalty function
%       numos   - number of ordered subsets
%       img0    - initial image
%       img1    - finished image
% output:
%       x       - reconstruction result
%       phis    - cost function values
%       rmsds   - root means squared differences with previous iteration or converged solution
%
% Based "Statistical image reconstruction for polyenergetic X-ray computed
% tomography" by Elbakri, Idris A and Fessler, Jeffrey A
%
% 2012-13 Meng Wu at Stanford University

if nargin < 10
    updateRate = 0.02;
end


if nargin < 11
    noFrames = 10;
end


figKey = round( 100 * rand());


fprintf('Generate contiueous reconstruction path for PWLS method with updata rate %0.3f%% ... \n', updateRate);
tRecon = tic;

% load operators for projection and regularization
[A, At, Aos, Atos, Os ] = loadPojectors( geom, numos );
[R, S, T ]  = loadPenaltyOperator( pfun, delta );

if R(img0) > R(img1)
    forwardPathSeeking = 1;
    fprintf('\tUsing forward mode. \n');
else
    forwardPathSeeking = 0;
    fprintf('\tUsing backward mode. \n');
end

x = img0;
diff = img1 - x;
validPixelsAll          = abs( diff ) >  dv;
totalNumUpdatePixels    = sum( validPixelsAll(:) );
maxNumUpdatePixels      = totalNumUpdatePixels * updateRate;
totalPathLength        = mean( abs(  diff(validPixelsAll(:)) ) );
noIterations            = round( totalPathLength / dv / updateRate / numos);
saveThreshold           = totalPathLength / ( noFrames + 1);


xs = zeros( [ size(x,1) size(x, 2) noFrames], 'single' );
phis = zeros(noFrames, 1);
regs = zeros(noFrames, 1);

fprintf('\tTotal %3i frames with step size %f, with about %i iterations.', noFrames, dv, noIterations );
hdg1 = '   itn       x(0)        ||x-x0||_1    ||x-x1||_1      PHI(x)         R(x)';
fprintf('\n%s'      , hdg1  );

xold    = x;
iframe  = 1;
itn     = 1;
totalAbsDistance = totalPathLength;
% xs(:,:, 1) = img0(:,:,ceil(end/2));

while totalAbsDistance >= saveThreshold && iframe <= noFrames
    
    %meanAbsErrorOld = totalAbsDistance;
    
    for isub = 1 : numos
        
        d = Aos( x, isub )- Os(y,isub) ;
        
        % gradient of data fidelity and regularization
        gradientWLS = Atos( Os(w,isub).* d, isub ) ;
        gradientPenalty = S(x);
        
        
        % normalize the gradiants
        gradientWLS = gradientWLS / std(gradientWLS(:));
        gradientPenalty = gradientPenalty / std(gradientPenalty(:));
        
        pixelSameDirection = gradientWLS .*  gradientPenalty > 0 ;
        
        
        if forwardPathSeeking
            lambda = gradientPenalty ./ ( abs( gradientWLS ) + 1e-3 );
            lambda( pixelSameDirection ) = lambda( pixelSameDirection ) * 100;
            % make sure the update pixel are in the right direction
            valid = abs( diff ) >  dv;
            valid = valid & ( ( lambda .* diff < 0 ) | pixelSameDirection);
        else
            lambda = gradientWLS ./ ( abs( gradientPenalty ) + 1e-3 );
            lambda( pixelSameDirection ) = lambda( pixelSameDirection ) * 100;
            % make sure the update pixel are in the right direction
            valid = abs( diff ) >  dv;
            valid = valid & ( lambda .* diff < 0 );
        end
        
        
        % rise up the percentage a little bit if only a few valid pixels
        q = min( maxNumUpdatePixels / sum(valid(:)), 0.5 );
        
        % find the 1-q quantile on the scores
        t = quantile( abs( lambda(valid(:)) ), 1 - q );
        
        % calculate the new differences
        diff = img1 - x;
        
        updatePixelMap = ( abs( lambda ) > t ) & valid;
        pixelUpdateFraction      = sum( updatePixelMap(:) ) / totalNumUpdatePixels;
        
%         if pixelUpdateFraction < updateRate / 5;
%             %updatePixelMap = abs( diff ) >  dv;
%             x = x + dv /  mean( abs( diff(:) ) ) * diff ;
%         else
%             x( updatePixelMap ) = x( updatePixelMap ) + dv * sign( diff( updatePixelMap ) );
%             %x( updatePixelMap ) = x( updatePixelMap ) - dv * sign( lambda( updatePixelMap ) );
%         end
        
        x( updatePixelMap ) = x( updatePixelMap ) + dv * sign( diff( updatePixelMap ) );
        x = x + dv * ( updateRate -  pixelUpdateFraction )  / mean( abs( diff(:) ) ) * diff ;
        
        
        
        
        totalAbsDistance    = mean( abs(  x(validPixelsAll(:)) - img1(validPixelsAll(:)) ) ) ;
        recentAbsUpdate     = mean( abs(  x(validPixelsAll(:)) - xold(validPixelsAll(:)) ) ) ;
        totalAbsUpdate      = mean( abs(  x(validPixelsAll(:)) - img0(validPixelsAll(:)) ) ) ;
        updatePercentage    = totalAbsUpdate / totalPathLength * 100 ;
        pixelUpdateFraction = sum( updatePixelMap(:) ) / totalNumUpdatePixels;
        
        % record a frame
        if  recentAbsUpdate >  saveThreshold
            
            gradientWLS = Atos( Os(w,isub).* d, isub ) ;
            gradientPenalty = S(x);
            
            gradientRatio =  - gradientWLS ./ gradientPenalty;
            valid = validPixelsAll & ~isnan( gradientRatio ) & ~isinf( gradientRatio ) & abs(gradientPenalty) > 1e-6;
            beta =  median( gradientRatio(valid(:)) );
            fprintf( '\n\tSave frame %i at beta = %f with percentage %f', iframe , beta, updatePercentage );
            
            xs(:,:,iframe ) = x(:,:,ceil(end/2));
            iframe = iframe + 1;
            xold = x;
        end
        
        if isub == 1
            phi = 0;
        end
        
        temp = Os(w, isub).* d.^2;
        phi = phi + sum( temp(:) ) / 2;
        
    end
    
    prnt = 0;
    if itn   <= 10       , prnt = 1; end
    if rem(itn, 10 ) == 0  , prnt = 1; end
    
    if prnt
        fprintf('\n%6g %13.3e'  , itn  , x(end/2,end/2,ceil(end/2)));
        fprintf(' %13.3e %13.3e', totalAbsUpdate  , totalAbsDistance );
        fprintf(' %13.3e %13.3e', phi  ,  R(x))
    end
    itn = itn + 1;
    
    if 0
        figure(figKey)
        subplot(231), imdisp( x(:,:,ceil(end/2)), [0.15 0.25] ); title( ['x at frame ' num2str(iframe - 1) ]);
        subplot(232), imdisp( lambda(:,:,ceil(end/2)), [-1000 1000] );  title( ['lambda at itn ' num2str(itn) ]);
        subplot(233), imdisp( updatePixelMap(:,:,ceil(end/2)) ); title( ['update map ' num2str(pixelUpdateFraction)]);
        subplot(234), imdisp( (x(:,:,ceil(end/2)) - img0(:,:,ceil(end/2))), [-0.02 0.02]); title( ['total update ' num2str(totalAbsUpdate)] );
        subplot(235), imdisp( (img1(:,:,ceil(end/2)) - x(:,:,ceil(end/2))), [-0.02 0.02] ); title( ['image error ' num2str(totalAbsDistance)] );
        subplot(236), imdisp( (x(:,:,ceil(end/2)) - xold(:,:,ceil(end/2))), [-dv dv] * 5 ); title( ['recent update ' num2str(recentAbsUpdate)] );
        pause(0.01)
    end
    
end

for jframe = iframe : noFrames
    
    wl = ( noFrames - jframe ) / ( noFrames - iframe + 1);
    wr = 1 - wl;
    xs(:,:,jframe) =  wl * xs(:,:,iframe - 1) + wr * img1(:,:,ceil(end/2));
    
end


tRecon = toc(tRecon);
fprintf('\nDone in %dmin %0ds.\n\n\n', floor(tRecon/60), round(mod(tRecon, 60)));



end

