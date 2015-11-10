function [xs, phis, regs] = reconPwlsPathSeeking2( y, w, geom, pfun, delta, numos, img0, img1, beta0, beta1, dv, p, noFrames )
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


if nargin < 14
    noFrames = 10;
end


figKey = round( 100 * rand());


fprintf('Generate contiueous reconstruction path for PWLS method with updata rate %0.3f%% ... \n', p);
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
maxNumUpdatePixels      = totalNumUpdatePixels * p;
pertangeUpdate          = p;
meanAbsError            = mean( abs(  diff(validPixelsAll(:)) ) ) / dv;
noIterations            = round( meanAbsError / p / ( numos + 1) );
saveThreshold           = meanAbsError / ( noFrames + 1);

a = A( ones( size(x), 'single' ) );
precom = At( w .* a );

xs = zeros( [ size(x)  noFrames], 'single' );
phis = zeros(noFrames, 1);
regs = zeros(noFrames, 1);

fprintf('\tTotal %3i frames with step size %f, with about %i iterations.', noFrames, dv, noIterations );
hdg1 = '   itn       x(0)          SD(x-x0)      SD(x-x1)      PHI(x)         R(x)';
fprintf('\n%s'      , hdg1  );

xold = x;
iframe = 0;
itn = 1;

while meanAbsError >= saveThreshold && iframe <= noFrames
    
    %meanAbsErrorOld = meanAbsError;
    
    for isub = 1 : numos
        
        d = Aos( x, isub )- Os(y,isub) ;
        
        % gradient of data fidelity and regularization
        gradientWLS = Atos( Os(w,isub).* d, isub ) ;
        gradientWLS = gradientWLS / std(gradientWLS(:));
        
        
        gradientPenalty = S(x);
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
        
        if pertangeUpdate < p / 10;
            updatePixelMap = abs( diff ) >  dv;
        end
        
        x( updatePixelMap ) = x( updatePixelMap ) + dv * sign( diff( updatePixelMap ) );
        %x( updatePixelMap ) = x( updatePixelMap ) - dv * sign( lambda( updatePixelMap ) );
        
        meanAbsError        = mean( abs(  x(validPixelsAll(:)) - img1(validPixelsAll(:)) ) ) / dv;
        recentAbsUpdate     = mean( abs(  x(validPixelsAll(:)) - xold(validPixelsAll(:)) ) ) / dv;
        totalAbsUpdate      = mean( abs(  x(validPixelsAll(:)) - img0(validPixelsAll(:)) ) ) / dv;
        pertangeUpdate      = sum( updatePixelMap(:) ) / totalNumUpdatePixels;
        
        % record a frame
        if  recentAbsUpdate >  saveThreshold
            
            d0 = mean((x(validPixelsAll(:)) - img0(validPixelsAll(:))).^2);
            
            d1 = mean((x(validPixelsAll(:)) - img1(validPixelsAll(:))).^2);
                        
            beta = (d1 * beta0 + d0 * beta1) / ( d0 + d1 );

            for jsub = 1 : numos
                
                d = Aos( x, jsub )- Os(y,jsub) ;
                
                gradient = Atos( Os(w,jsub).* d, jsub ) + beta * S(x);
                
                dx =  gradient ./ ( precom + beta * T(x) ) ;
                
                x = x - dx;
                
                x( x < 0 ) = 0;
                x( isnan(x) ) = 0;
                x( isinf(x) ) = 0;
            end
            
            
            xs(:,:,iframe + 1) = x;
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
        fprintf(' %13.3e %13.3e', std(x(validPixelsAll(:)) - img0(validPixelsAll(:)))  , std(x(validPixelsAll(:)) - img1(validPixelsAll(:))));
        fprintf(' %13.3e %13.3e', phi  ,  R(x))
    end
    itn = itn + 1;
    
    
    if 0
        figure(figKey)
        subplot(231), imdisp( x', [0.15 0.25] ); title( ['x at frame ' num2str(iframe) ]);
        subplot(232), imdisp( lambda', [-1000 1000] );  title( ['lambda at itn ' num2str(itn) ]);
        subplot(233), imdisp( single( updatePixelMap' ) ); title( ['update map ' num2str(pertangeUpdate)]);
        subplot(234), imdisp( (x - img0)', [-0.02 0.02]); title( ['total update ' num2str(totalAbsUpdate)] );
        subplot(235), imdisp( (img1 - x)', [-0.02 0.02] ); title( ['image error ' num2str(meanAbsError)] );
        subplot(236), imdisp( (x - xold)', [-dv dv] * 5 ); title( ['recent update ' num2str(recentAbsUpdate)] );
        pause(0.01)
    end
    
    
    
end

tRecon = toc(tRecon);
fprintf('\nDone in %dmin %0ds.\n\n\n', floor(tRecon/60), round(mod(tRecon, 60)));



end

