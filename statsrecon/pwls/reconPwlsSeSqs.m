function [x, phis, rmsds] = reconPwlsSeSqs( y, w, geom, beta, pfun, itnlim, delta, numos, img0, imgOpt  )
% Penalized-Weight-Least-Squares recosntruction for single energy sinograms
% using Nonuniform SQS algorithm
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
% Based "Statistical image reconstruction for polyenergetic X-ray computed
% tomography" by Elbakri, Idris A and Fessler, Jeffrey A
%
% 2012-13 Meng Wu at Stanford University

fprintf('Reconstructing with keV data using PWLS SQS... \n');
tRecon = tic;

% use FBP to compute initial image
if nargin < 9
    img0 = reconFBP( y, geom );
end

stopCrt         = 1e-5;

% load operators for projection and regularization
[A, At, Aos, Atos, Os ] = loadPojectors( geom, numos );
[R, S, T ]  = loadPenaltyOperator( pfun, delta );

x = img0;
a = A( ones( size(x), 'single' ) );
precom = At( w .* a );

fprintf('\nbeta  = %11.2e', beta );
fprintf('\titnlim = %10g', itnlim);
hdg1 = '   itn       x(0)          PHI(x)        b*R(x)        RMSD';
fprintf('\n%s'      , hdg1  );

phis = zeros(itnlim, 1);
rmsds = zeros(itnlim, 1);

for itn = 1 : itnlim
    
    xold = x;
    
    for isub = 1 : numos
        
        d = Aos( x, isub )- Os(y,isub) ;
        
        numerator = Atos( Os(w,isub).* d, isub ) + beta * S(x);
        
        denominator = precom + beta * T(x);
        
        x = x - numerator ./ denominator;
        x( x < 0 ) = 0;
        x( isnan(x) ) = 0;
        x( isinf(x) ) = 0;
        
        
        if isub == 1
            phi = 0;
        end
        
        temp = Os(w, isub).* d.^2;
        phi = phi + sum( temp(:) ) / 2;
        
    end
    
    phi = phi + beta * R(x);
    
    if nargin < 10
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


tRecon = toc(tRecon);
fprintf('\nDone in %dmin %0ds.\n\n\n', floor(tRecon/60), round(mod(tRecon, 60)));




end

