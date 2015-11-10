function [xs, betas, phis, regs] = reconPwlsTruePathSeeking( y, w, geom, pfun, delta, psOS, optOS, ...
    img0, img1, beta0, beta1, dv0, updateRate, noFrames, numOptItn )
% Penalized-Weight-Least-Squares recosntruction for single energy sinograms
% using Nonuniform SQS algorithm
% input:
%       y       - log sinogram
%       w       - weights for pwls
%       geom    - system geometry
%       pfun    - roughness penalty function (default 'huber')
%       delta   - parameter for penalty function
%       psOS   - number of ordered subsets
%       img0    - initial image
%       img1    - finished image
% output:
%       x       - reconstruction result
%       phis    - cost function values
%       rmsds   - root means squared differences with previous iteration or converged solution
%
%
% 2012-13 Meng Wu at Stanford University

if nargin < 10
    updateRate = 0.02;
end


if nargin < 11
    noFrames = 10;
end

incrementRatio  = exp( log( beta1 / beta0 ) /  noFrames );

fprintf('Path seeking for PWLS method using ratio of gradients ... \n');
tRecon = tic;

k  = coneTruncatedSlices( geom );

% load operators for projection and regularization
[A, At, Aos, Atos, Os ] = loadPojectors( geom, psOS );
[~, ~, Aos2, Atos2, Os2 ] = loadPojectors( geom, optOS );

[R, S, T ]  = loadPenaltyOperator( pfun, delta );

a = A( ones( size(img0), 'single' ) );
precom = At( w .* a );
precom = extendVoi( precom, k );

if beta0 < beta1
    forwardPathSeeking = 1;
    fprintf('\tUsing forward mode. \n');
else
    forwardPathSeeking = 0;
    fprintf('\tUsing backward mode. \n');
end

x = img0;
diff = img1 - x;
diff = extendVoi( diff, k + 1 );

validPixelsAll          = abs( diff ) >  dv0;
totalNumUpdatePixels    = sum( validPixelsAll(:) );
totalPathLength         = mean( abs( diff(:) ) );

% determine the step size
dv                      = mean( abs(  diff(validPixelsAll(:)) ) ) / ( updateRate * noFrames * psOS );
noSubIterations         = ceil( dv / dv0 );
dv                      = dv / noSubIterations;

% determin the number of iteration
noIterations            = noFrames * noSubIterations;

xs = zeros( [ size(x,1) size(x, 2) noFrames], 'single' );
betas = zeros( 1, noFrames );
phis = zeros(noFrames, 1);
regs = zeros(noFrames, 1);

fprintf('\tThe updata rate is %0.3f%% with step size %f.\n', updateRate, dv);
fprintf('\tTotal %3i frames need %i iterations.\n', noFrames, noIterations * ( 1 + numOptItn ) );
fprintf('\tAdditiona %i sub-iteration(s) to improve accuracy. \n', numOptItn )
fprintf('   itn       x(0)         ||x-x0||_1    ||x-x1||_1      PHI(x)         R(x)' );

iframe  = 1;

beta = beta0;
betaIncrement = beta0;
% beta = estimateBetaKKT( At( w .* ( A(x) - y)), S(x), validPixelsAll, 100, forwardPathSeeking );

for itn = 1 : noIterations
    
    gradientWLSAcc = 0;
    for isub = 1 : psOS
        
        
        % gradient of data fidelity and regularization
        d = Aos( x, isub )- Os(y,isub);
        gradientWLS         = Atos( Os(w,isub).* d, isub ) ;
        gradientPenalty 	= S(x);
        
        gradientWLSAcc = gradientWLSAcc + gradientWLS;
        
        % a gradient descent step
        dx =  (gradientWLS + beta * gradientPenalty ) ./ ( precom + beta * T(x) ) ;
        
        % calculate the new differences
        diff = img1 - x;
        % calculate the update the ratio of gradients
        dx = dx + gradientRatioUpdata( gradientWLS, gradientPenalty, diff, ...
            forwardPathSeeking, dv, updateRate, totalNumUpdatePixels );
        
        % make the update here
        x = x - dx;
        x( x < 0 ) = 0;
        x( isnan(x) ) = 0;
        x( isinf(x) ) = 0;
        
        if isub == 1
            phi = 0;
        end
        
        temp = Os(w, isub).* d.^2;
        phi = phi + sum( temp(:) ) / 2;
        
    end
    
    % using the KKT condition to estimate the tuning parameter value
    beta = estimateBetaKKT( gradientWLSAcc, gradientPenalty, validPixelsAll, beta, forwardPathSeeking );
    beta = sqrt( beta * betaIncrement );
    
    for iopt = 1 : numOptItn
        for isub = 1 : optOS
            numerator = Atos2( Os2(w,isub).* (Aos2( x, isub )- Os2(y,isub)), isub ) + beta * S(x);
            x = x - numerator ./ ( precom + beta * T(x) );
            x( x < 0 ) = 0;
            x( isnan(x) ) = 0;
            x( isinf(x) ) = 0;
            % x = extendVoi( x, k );
        end
    end

    totalAbsDistance    = mean( abs(  x(validPixelsAll(:)) - img1(validPixelsAll(:)) ) ) ;
    totalAbsUpdate      = mean( abs(  x(validPixelsAll(:)) - img0(validPixelsAll(:)) ) ) ;
    updatePercentage    = totalAbsUpdate / totalPathLength * 100 ;
    
    fprintf('\n%6g %13.3e'  , itn  , x(round(end/2), round(end/2), ceil(end/2)));
    fprintf(' %13.3e %13.3e', totalAbsUpdate  , totalAbsDistance );
    fprintf(' %13.3e %13.3e', phi  ,  R(x))
    
    % record a frame
    if  mod( itn , noSubIterations ) == 0
        fprintf( '\tSave frame %i with percentage %.2f, beta = %.3g.', iframe , updatePercentage,  beta );
        xs(:,:,iframe) = x(:,:,ceil(end/2));
        betas( iframe ) = beta;
        phis( iframe ) = phi;
        regs( iframe ) = R(x);
        
        iframe = iframe + 1;
        
        betaIncrement = betaIncrement * incrementRatio;
        
    end
    
end

tRecon = toc(tRecon);
fprintf('\nDone in %dmin %0ds.\n\n\n', floor(tRecon/60), round(mod(tRecon, 60)));

end


