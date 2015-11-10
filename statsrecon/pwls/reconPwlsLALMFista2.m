function [x, phis, rmsds] = reconPwlsLALMFista2( y, w, geom, beta, pfun, itnlim, delta, numos, imgOpt  )
% Penalized-Weight-Least-Squares recosntruction for single energy sinograms
% using Nonuniform SQS algorithm
% input:
%       y       - log sinogram
%       w       - weights for pwls
%       geom    - system geometry
%       beta    - weight on penalty function
%       pfun    - roughness penalty function (default 'huber')
%       itnlim    - maximum number of iterations
%       isFlat  - is flat penal
% output:
%       x       - reconstruction result
%       phis    - cost function values
%
% Based Nien, Hung, and Jeffrey A. Fessler. “Fast X-Ray CT Image
% Reconstruction Using the Linearized Augmented Lagrangian Method with
% Ordered Subsets.” Optimization and Control; Learning; Machine Learning.
% arXiv Preprint arXiv:1402.4381 (February 18, 2014): 21. 
%
% Eqn. (33) 
%
% 2012-13 Meng Wu at Stanford University

fprintf('Reconstructing with single energy data using PWLS LALM ... \n');
tRecon = tic;


if nargin < 4
    beta = 1e3;
end

if nargin < 5
    pfun = 'huber';
end

if nargin < 6
    itnlim = 200;
end

if nargin < 7
    delta           = 0.001;
end

if nargin < 8
    numos           = 24;
end

stopCrt         = 1e-5;
showCostFunc    = false;
k               = 1;

rho = 0.1;

% select projection geometry model
Aos     = @( x, is )forwardProjectMex( x, geom, numos, is, 'proj,dd' );
Atos    = @( x, is )backProjectMex( x, geom, numos, is, 'back,dd' );
A       = @( x )forwardProjectMex( x, geom, 1, 0, 'proj,dd' );
At      = @( x )backProjectMex( x, geom, 1, 0, 'back,dd' );
Os      = @( x, is)orderedSubset(x, numos , is);

if strcmpi(pfun, 'huber')
    R       = @( x )huberPenalty( x, 0, delta );
    S       = @( x )huberPenalty( x, 1, delta );
    T       = @( x )huberPenalty( x, 2, delta );
elseif strcmpi(pfun, 'quad')
    R       = @( x )quadPenalty( x, 0 );
    S       = @( x )quadPenalty( x, 1 );
    T       = @( x )quadPenalty( x, 2 );
elseif strcmpi(pfun, 'hyper')
    R       = @( x )hyperbolaPenalty( x, 0, delta );
    S       = @( x )hyperbolaPenalty( x, 1, delta );
    T       = @( x )hyperbolaPenalty( x, 2, delta );
elseif strcmpi(pfun, 'aniso')
    R       = @( x )anisotropicPenalty( x, 0, delta );
    S       = @( x )anisotropicPenalty( x, 1, delta );
    T       = @( x )anisotropicPenalty( x, 2, delta );
else
    error('unknow penalty function');
end

% use FBP to compute initial image
x = reconFBP( y, geom );
x = extendVoi( x, k );

a = A( ones( size(x), 'single' ) );
precom = At( w .* a );
t = 1 / max( precom(:) );

fprintf('\nbeta  = %11.2e, rho  = %11.2e\n', beta, rho );
fprintf('\titnlim = %10g', itnlim);
hdg1 = '   itn       x(0)          PHI(x)        b*R(x)        RMSD';
fprintf('\n%s'      , hdg1  );

phis = zeros(itnlim, 1);
rmsds = zeros(itnlim, 1);

g = zeros( size(x), 'single');
for itn = 1 : itnlim
    
    xold = x;

    % step 3: d^{(k+1)} = d^{(k)} - A x^{k+1} + u^{(k+1)} 
    r = A( x ) - y;
    dL = At( w.* r  )  ;
    g = rho / (rho + 1) * dL + 1/( rho + 1 ) * g;
    
    % step 1: s^{(k+1 )} = rho * dL( x^{(k)} ) + ( 1 - rho) g^{(k)} 
    s = rho * dL + ( 1 - rho ) * g ;
    s( isnan(s) ) = 0;
    s( isinf(s) ) = 0;
    
    % step 2: x^{(k+1)} = prox_{ rho^-1 t }( x^{k} - ( rho^-1 t ) s^{k+1})
    % TODO: implemented by FISTA
    
    for j = 1:10
        alpha = 1 / rho * t ;
        
        x = ( x + alpha * beta * S(x) - alpha * s  ) ./ ( alpha * beta *  T(x)  + 1 ); 
        x( x < 0 ) = 0;
        x( isnan(x) ) = 0;
        x( isinf(x) ) = 0;
    end
    
    
    phi = sum(  w(:).* r(:).^2 ) / 2 + beta * R(x);
    if nargin < 9
        rmsd = sqrt( mean( (x(:) - xold(:)).^2 ) );
    else
        if ndims( x ) == 3
            a = x(:,:,end/2);
            b = imgOpt(:,:,end/2);
        else
            a = x;
            b = imgOpt;
        end
        rmsd = sqrt( mean( (a(:) - b(:)).^2 ) );
    end
    phis(itn) = phi;
    rmsds(itn) = rmsd;
    
    prnt = 0;
    if itn   <= 10       , prnt = 1; end
    if itn   >= itnlim-10, prnt = 1; end
    if rem(itn,10) == 0  , prnt = 1; end
    
    if prnt
        fprintf('\n%6g %13.3e'  , itn  , x(end/2,end/2,ceil(end/2)));
        fprintf(' %13.3e %13.3e %13.3e', phi  , beta * R(x), rmsd)
    end
    
    if itn > 5 && rmsd < stopCrt
        break;
    end
    
end


if showCostFunc
    figure('Name', 'cost function values', 'NumberTitle', 'off');
    semilogy( 1:itn, phis(1:itn) );
    ylabel '\phi(x)'; xlabel 'No. of iteration';
    
    figure('Name', 'RMSD values', 'NumberTitle', 'off');
    semilogy( 1:itn, rmsds(1:itn) );
    ylabel 'RMDS'; xlabel 'No. of iteration';
    
end

tRecon = toc(tRecon);
fprintf('\nDone in %dmin %0ds.\n', floor(tRecon/60), round(mod(tRecon, 60)));



end

