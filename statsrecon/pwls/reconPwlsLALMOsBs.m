function [x, phis, rmsds] = reconPwlsLALMOsBs( y, w, geom, beta, pfun, itnlim, delta, numos, imgOpt  )
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

fprintf('Reconstructing with single energy data using PWLS LALM OS BS... \n');
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

numas       =  4;
rho         = 0.1;


maps = chessboardMaps( geom, numas  );

% load operators for projection and regularization
[A, At, Aos, Atos, Os ] = loadPojectors( geom, numos );
[R, S, T ]  = loadPenaltyOperator( pfun, delta );


% use FBP to compute initial image
x = reconFBP( y, geom );
x = extendVoi( x, k );

a = A( ones( size(x), 'single' ) );
precom = At( w .* a );

fprintf('\nbeta  = %11.2e, rho  = %11.2e\n', beta, rho );
fprintf('\titnlim = %10g', itnlim);
hdg1 = '   itn       x(0)          PHI(x)        b*R(x)        RMSD';
fprintf('\n%s'      , hdg1  );

phis = zeros(itnlim, 1);
rmsds = zeros(itnlim, 1);

g = zeros( size(x), 'single');

u = cell(numas, 1);
v = cell(numas, 1);


for itn = 1 : itnlim
    
    for isub = 1 : numos
        
        
        xold = x;
        
        % step 3: d^{(k+1)} = d^{(k)} - A x^{k+1} + u^{(k+1)} 
        parfor jsub = 1 : numas 
            u{jsub} =  Aos( x, isub, maps{jsub} ) ;
        end
        
        r = - Os( y, isub );
        for jsub = 1:numas
           r = r + u{jsub}; 
        end
        

        parfor jsub = 1 : numas 
            v{jsub} = Atos( Os( w, isub) .* r, isub, maps{jsub} ) ;
        end
        
        dL = zeros( size(x), 'single');
        for jsub = 1 : numas 
            dL = dL + v{jsub} ;
        end

        g = rho / (rho + 1) * dL + 1/( rho + 1 ) * g;
        
        % step 1: s^{(k+1 )} = rho * dL( x^{(k)} ) + ( 1 - rho) g^{(k)}
        s = rho * dL + ( 1 - rho ) * g ;
        s( isnan(s) ) = 0;
        s( isinf(s) ) = 0;
        
        % step 2: x^{(k+1)} = prox_{ rho^-1 t }( x^{k} - ( rho^-1 t ) s^{k+1})
        % TODO: implemented by FISTA
        
        for j = 1:10
            x = x - ( beta * S(x) + rho * precom .* ( x - xold ) + s  ) ./ ( beta * T(x) + rho * precom );
            x( x < 0) = 0;
            x( isnan(x) ) = 0;
            x( isinf(x) ) = 0;
        end
        
        if isub == 1
            phi = 0;
        end
        
        
        temp = Os(w, isub).* r.^2;
        phi = phi + sum( temp(:) ) / 2;
        
    end
    
    phi = phi + beta * R(x);
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

