function [x, phis, rmsds] = reconPwlsSeNesterovSqsResample( y, w, geom, beta, pfun, itnlim, delta, numos, img0, imgOpt )
% Penalized-Weight-Least-Squares recosntruction for single energy sinograms
% using Nestrov's 05 + Noneuniform SQS algorithm
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
%Kim, Donghwan, Sathish Ramani, and Jeffrey A. Fessler. “Accelerating X-Ray
%CT Ordered Subsets Image Reconstruction with Nesterov’s FIrst-Order
%Method.” In The 12th International Meeting on Fully Three-Dimensional
%Image Reconstruction in Radiology and Nuclear Medicine, 22 – 25, 2013.
%
% 2013 - 2014  Meng Wu at Stanford University

fprintf('Reconstructing with keV data using PWLS Nesterov05 SQS ... ');
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
    delta           = 1;
end

if nargin < 8
    numos           = 41;
end

% use FBP to compute initial image
if nargin < 9
    img0 = reconFBP( y, geom );
end

stopCrt         = 1e-10;
nrestart        = 20;
% load operators for projection and regularization
[A, At, Aos, Atos, Os ] = loadPojectors( geom, numos );
[R, S, T ]  = loadPenaltyOperator( pfun, delta );

x = img0;
a = A( ones( size(x), 'single' ) );
precom = At( w .* a );

w = w .* single( rand( size(w) ) > (1 - delta) );

fprintf('\nbeta  = %11.2e', beta );
fprintf('\titnlim = %10g', itnlim);
hdg1 = '   itn       x(0)          PHI(x)        b*R(x)        RMSD';

fprintf('\n%s'      , hdg1  );

phis = zeros(itnlim, 1);
rmsds = zeros(itnlim, 1);

%initialized for Nesterov's algo
t = 1;
z = x;
s = x;

for itn = 1 : itnlim
    
    xold = x;
    
    for isub = 1 : numos
        
        d = Aos( z, isub )- Os(y,isub) ;
        
        gradient = Atos( Os(w,isub).* d, isub ) + beta * S(z);
        
        du =  gradient ./ ( precom + beta * T(x) ) ;
        
        x = z - du;
        s = s - t * du;
        
        t = ( 1 + sqrt( 1 + 4 * t^2 ) ) / 2;
        z = ( 1 - 1/t ) * x + (1/t) * s;
        
        x( x < 0 ) = 0;
        x( isnan(x) ) = 0;
        x( isinf(x) ) = 0;
        
        s( s < 0 ) = 0;
        s( isnan(s) ) = 0;
        s( isinf(s) ) = 0;
        
        z( z < 0 ) = 0;
        z( isnan(z) ) = 0;
        z( isinf(z) ) = 0;
        
        
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
    
    if mod( itn, nrestart ) == 0 % restart nestrov
        t = 1;
        z = x;
        s = x;
    end
    
end

tRecon = toc(tRecon);
fprintf('\nDone in %dmin %0ds.\n\n\n', floor(tRecon/60), round(mod(tRecon, 60)));


end
