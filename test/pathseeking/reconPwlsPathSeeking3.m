function [xs, phis, regs] = reconPwlsPathSeeking3( y, w, geom, pfun, delta, numos, img0, img1, dv, p )
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
    p = 0.02;
end

figKey = round( 100 * rand());


fprintf('Generate contiueous reconstruction path for PWLS method with updata rate %0.3f%% ... ', p);
tRecon = tic;

% load operators for projection and regularization
[A, At, Aos, Atos, Os ] = loadPojectors( geom, numos );
[R, S, T ]  = loadPenaltyOperator( pfun, delta );




if R(img0) > R(img1)
    forwardPathSeeking = 1;
else
    forwardPathSeeking = 0;
end

x = img0;
diff = img1 - x;
validPixelsAll          = abs( diff ) >  dv;
totalNumUpdatePixels    = sum( validPixelsAll(:) );
maxNumUpdatePixels      = totalNumUpdatePixels * p;
meanAbsError            = mean( abs(  diff(validPixelsAll(:)) ) ) / dv;
noFrames                = round( meanAbsError + 1);
noIterations            = round( noFrames / p / numos );
meanAbsErrorOld = meanAbsError + 1;


xs = zeros( [ size(x)  noFrames], 'single' );
phis = zeros(noFrames, 1);
regs = zeros(noFrames, 1);

fprintf('\n\tTotal %3i frames with step size %f, with about %i iterations.', noFrames, dv, noIterations );
hdg1 = '   itn       x(0)          SD(x-x0)      SD(x-x1)      PHI(x)         R(x)';
fprintf('\n%s'      , hdg1  );

xold = x;
iframe = 1;
itn = 1;
updatePixelMap = validPixelsAll;

d = A( x )- y ;
% gradient of data fidelity and regularization
gradientWLS = At( w .* d ) ;
gradientWLS( abs( gradientWLS )  < 1e-8  ) = 1e-8;
gradientWLS = gradientWLS / std(gradientWLS(:));


while iframe < noFrames
    
    meanAbsErrorOld = meanAbsError;
    
    gradientPenalty = S(x);
    gradientPenalty( abs( gradientPenalty )  < 1e-8  ) = 1e-8;
    gradientPenalty = gradientPenalty / std(gradientPenalty(:));
    
    
    
    if forwardPathSeeking
        lambda = gradientPenalty ./ abs( gradientWLS  );
    else
        lambda = gradientWLS ./ abs( gradientPenalty );
    end
    
    if sum( updatePixelMap(:) ) < maxNumUpdatePixels * 0.1
        lambda = lambda - diff / dv;
    end
    
    
    lambda( ~validPixelsAll ) = 0;
    %lambda( lambda .* diff > 0 ) = 0;
    
    % make sure the update pixel are in the right direction
    %valid = abs( diff ) >  dv;
    %valid = valid & ( lambda .* diff < 0 );
    
    
    % rise up the percentage a little bit if only a few valid pixels
    %q = min( maxNumUpdatePixels / sum(valid(:)), 0.5 );
    
    % find the 1-q quantile on the scores
    t = quantile( abs( lambda(validPixelsAll(:)) ), 1 - p );
    
    % get a window for displaying
    if iframe == 1
        if forwardPathSeeking
            t0 = t * 3;
        else
            t0 = t;
        end
    end
    
    updatePixelMap = ( abs( lambda ) > t );
    x( updatePixelMap ) = x( updatePixelMap ) + dv * sign( diff( updatePixelMap ) );
    % x( updatePixelMap ) = x( updatePixelMap ) - dv * sign( lambda( updatePixelMap ) );
    
    % calculate the new differences
    diff = img1 - x;
    
    meanAbsError        = mean( abs(  x(validPixelsAll(:)) - img1(validPixelsAll(:)) ) ) / dv;
    recentAbsUpdate     = mean( abs(  x(validPixelsAll(:)) - xold(validPixelsAll(:)) ) ) / dv;
    totalAbsUpdate      = mean( abs(  x(validPixelsAll(:)) - img0(validPixelsAll(:)) ) ) / dv;
    pertangeUpdate      = sum( updatePixelMap(:) ) / totalNumUpdatePixels;
    
    if  recentAbsUpdate >  1
        
        d = A( x )- y ;
        % gradient of data fidelity and regularization
        gradientWLS = At( w .* d ) ;
        gradientWLS( abs( gradientWLS )  < 1e-8  ) = 1e-8;
        gradientWLS = gradientWLS / std(gradientWLS(:));
        
        
        xs(:,:,iframe) = x;
        iframe = iframe + 1;
        xold = x;
        
    end
    
    
    temp = w.* d.^2;
    phi = sum( temp(:) ) / 2;
    
    if 1
        figure(figKey)
        subplot(231), imdisp( x', [0.15 0.25] ); title( ['x at frame ' num2str(iframe) ]);
        subplot(232), imdisp( lambda', [-t0 t0] );  title( ['lambda at itn ' num2str(itn) ]);
        subplot(233), imdisp( single( updatePixelMap' ) ); title( ['update map ' num2str(pertangeUpdate)]);
        subplot(234), imdisp( (x - img0)', [-0.02 0.02]); title( ['total update ' num2str(totalAbsUpdate)] );
        subplot(235), imdisp( (img1 - x)', [-0.02 0.02] ); title( ['image error ' num2str(meanAbsError)] );
        subplot(236), imdisp( (x - xold)', [-dv dv] * 5 ); title( ['recent update ' num2str(recentAbsUpdate)] );
        pause(0.01)
    end
    
    
    
    prnt = 0;
    if itn   <= 10       , prnt = 1; end
    if rem(itn, 10 ) == 0  , prnt = 1; end
    
    if prnt
        fprintf('\n%6g %13.3e'  , itn  , x(end/2,end/2,ceil(end/2)));
        fprintf(' %13.3e %13.3e', std(x(validPixelsAll(:)) - img1(validPixelsAll(:)))  , std(x(validPixelsAll(:)) - img0(validPixelsAll(:))));
        fprintf(' %13.3e %13.3e', phi  ,  R(x))
    end
    itn = itn + 1;
    
end

xs = xs(:,:,1:iframe-1);
phis = phis(1:iframe-1);
regs = regs(1:iframe-1);

tRecon = toc(tRecon);
fprintf('\nDone in %dmin %0ds.\n\n\n', floor(tRecon/60), round(mod(tRecon, 60)));



end

