function [xs, betas, phis, regs] = reconPwlsPathSeekingADMM( y, w, geom, pfun, delta, numos, img0, img1, ...
    beta0, beta1, noFrames, noSubIterations )
% Penalized-Weight-Least-Squares recosntruction for single energy sinograms
% using linearized augmented largragian + ordered subset
% input:
%       y       - log sinogram
%       w       - weights for pwls
%       geom    - system geometry
%       beta    - weight on penalty function
%       pfun    - roughness penalty function (default 'huber')
%       itnlim  - maximum number of iterations
%       delta   - parameter for penalty function
%       numos   - number of ordered subsets
%       img0    - initial image
%       imgOpt  - converged image
% output:
%       x       - reconstruction result
%       phis    - cost function values
%       rmsds   - root means squared differences with previous iteration or converged solution
%
% Based Nien, Hung, and Jeffrey A. Fessler. “Fast X-Ray CT Image
% Reconstruction Using the Linearized Augmented Lagrangian Method with
% Ordered Subsets.” Optimization and Control; Learning; Machine Learning.
% arXiv Preprint arXiv:1402.4381 (February 18, 2014): 21.
%
% Eqn. (33)
%
% 2012-13 Meng Wu at Stanford University

fprintf('Path seeking using for PWLS reonstruction using ADMM ... \n');
tRecon = tic;

% use FBP to compute initial image
if nargin < 9
    img0 = reconFBP( y, geom );
end

k  = coneTruncatedSlices( geom );

% load operators for projection and regularization
[A, At, Aos, Atos, Os ] = loadPojectors( geom, numos );
[R, S, T ]  = loadPenaltyOperator( pfun, delta );


a = A( ones( size(img0), 'single' ) );
precom = At( w .* a );
precom( precom < 1000 ) = 1000;
precom = extendVoi( precom, k );

beta    = beta0;
x       = img0;

% determin the number of iteration
totalPathLength  = mean( abs( img1(:) - img0(:) ));
noIterations     = noFrames * noSubIterations;
incrementRatio  = exp( log( beta1 / beta0 ) /  noFrames );

fprintf('\tTotal %3i frames need %i iterations, beta increment rate is %f.\n', noFrames, noIterations, incrementRatio );
if beta0 < beta1
    forwardPathSeeking = 1;
    fprintf('\tUse forward mode. ');
else
    forwardPathSeeking = 0;
    fprintf('\tUse backward mode. ');
end
fprintf('\tUse %i sub-iteration(s) to improve accuracy.\n', noSubIterations )
fprintf('   itn       x(0)         ||x-x0||_1    ||x-x1||_1      PHI(x)         R(x)' );

xs = zeros( [ size(x,1) size(x, 2) noFrames], 'single' );
phis = zeros(noFrames, 1);
regs = zeros(noFrames, 1);
betas = zeros(noFrames, 1);

iframe = 1;
rho = 0.05;

g = zeros( size(x), 'single');

% lastSaveInt = 0;
% while iframe <= noFrames

for itn = 1 : noIterations
    
    %g = zeros( size(x), 'single');
    for isub = 1 : numos
        
        xold = x;
        
        % step 3: d^{(k+1)} = d^{(k)} - A x^{k+1} + u^{(k+1)}
        r = Aos( x, isub ) - Os( y, isub );
        dL = Atos( Os( w, isub) .* r, isub  )  ;
        g = rho / (rho + 1) * dL + 1/( rho + 1 ) * g;
        
        % step 1: s^{(k+1 )} = rho * dL( x^{(k)} ) + ( 1 - rho) g^{(k)}
        s = rho * dL + ( 1 - rho ) * g ;
        s( isnan(s) ) = 0;
        s( isinf(s) ) = 0;
        
        % step 2: x^{(k+1)} = prox_{ rho^-1 t }( x^{k} - ( rho^-1 t ) s^{k+1})
        % TODO: implemented by FISTA
        l = S(x);
        for j = 1:5
            dx = ( beta * S(x) + rho * precom .* ( x - xold ) + s  ) ./ ( beta * T(x) + rho * precom );
            
            % direction of gradient type updates
            if mod( itn + isub , noSubIterations ) == 1
                if forwardPathSeeking 
                    dx( dx .* l < 0 ) = 0;
                else
                    dx( dx .* l > 0 ) = 0;
                end
            end
            x = x - dx;
            x( x < 0) = 0;
            x( isnan(x) ) = 0;
            x( isinf(x) ) = 0;
        end
        %x = extendVoi( x, k );
        
        if isub == 1
            phi = 0;
        end
        
        temp = Os(w, isub).* r.^2;
        phi = phi + sum( temp(:) ) / 2;
        
    end
    
    totalAbsDistance    = mean( abs(  x(:) - img1(:) ) ) ;
    totalAbsUpdate      = mean( abs(  x(:) - img0(:) ) ) ;
    updatePercentage    = totalAbsUpdate / totalPathLength * 100 ;
    
    fprintf('\n%6g %13.3e'  , itn  , x(round(end/2),round(end/2),ceil(end/2)));
    fprintf(' %13.3e %13.3e', totalAbsUpdate  , totalAbsDistance );
    fprintf(' %13.3e %13.3e', phi  ,  R(x))
    
    if  mod( itn , noSubIterations ) == 0
        fprintf( '\tSave frame %i with percentage %.2f, beta = %.3g.', iframe , updatePercentage,  beta );
        xs(:,:,iframe) = x(:,:,ceil(end/2));
        betas( iframe ) = beta;
        phis( iframe ) = phi;
        regs( iframe ) = R(x);
        
        iframe = iframe + 1;
        beta = beta * incrementRatio;
    end
    
    
end


tRecon = toc(tRecon);
fprintf('\nDone in %dmin %0ds.\n\n\n', floor(tRecon/60), round(mod(tRecon, 60)));



end

