function [x, phis, rmsds] = reconPwlsSeNesterovNusqs( y, w, geom, beta, pfun, itnlim, delta, numos, img0, imgOpt )
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
% Kim, Donghwan, Debashish Pal, Jean-Baptiste Thibault, and Jeffrey A.
% Fessler. “Accelerating Ordered Subsets Image Reconstruction for X-Ray CT
% Using Spatially Nonuniform Optimization Transfer.” IEEE Transactions on
% Medical Imaging 32, no. 11 (November 07, 2013): 1965–78. doi:10.1109/TMI.2013.2266898. 
%
% Kim, Donghwan, Sathish Ramani, and Jeffrey A. Fessler. “Accelerating X-Ray
% CT Ordered Subsets Image Reconstruction with Nesterov’s FIrst-Order
% Method.” In The 12th International Meeting on Fully Three-Dimensional
% Image Reconstruction in Radiology and Nuclear Medicine, 22 – 25, 2013.
%
% 2013 - 2014  Meng Wu at Stanford University

fprintf('Reconstructing with keV data using PWLS Nesterov05 NUSQS ... \n');
tRecon = tic;

% use FBP to compute initial image
if nargin < 9
    img0 = reconFBP( y, geom );
end

stopCrt         = 1e-5;
epsilon         = 1e-3;
nloop           = 3;
nfix            = 7;

% load operators for projection and regularization
[A, At, Aos, Atos, Os ] = loadPojectors( geom, numos );
[R, S, ~, H ]  = loadPenaltyOperator( pfun, delta );


x = img0;
% initalization for u0 using edge and intensity of FBP
v = zeros( size(x), 'single');
v = v + abs( imageFilter2D( x, [1 2 1; 0 0 0; -1 -2 1] ) );
v = v + abs( imageFilter2D( x, [1 0 -1; 2 0 -2; 1 0 -1] ) );
v = v * 2 + abs( x );
v( v < epsilon ) = epsilon;


dL =  ( At( w.* A(v) ) ) ./ v;

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
    
    if mod( itn, nloop ) == 1 && itn < nfix
        dLtilde = dL;
        dRtilde = beta *  H( v );
    elseif mod( itn, nloop ) == nloop - 1 && itn <= nfix - 2
        xref = x;
    elseif mod( itn, nloop ) == 0 && itn <= nfix - 1
        v = abs( x - xref );
        v( v < epsilon) = epsilon;
        v = dynamicRangeAdjustment( v, 40, 0.05  );
        dLtilde = ( At( w.* A(v) ) ) ./ v;
    end
    
    xold = x;
    
    for isub = 1 : numos
        
        d = Aos( z, isub )- Os(y,isub) ;
        gradient = Atos( Os(w,isub).* d, isub ) + beta * S(z);
        du =  gradient ./ ( dLtilde + dRtilde ) ;
        
        x = z - du;
        s = s - t * du;
        
        x( x < 0 ) = 0;
        x( isnan(x) ) = 0;
        x( isinf(x) ) = 0;
        

        s( s < 0 ) = 0;
        s( isnan(s) ) = 0;
        s( isinf(s) ) = 0;
        
        t = ( 1 + sqrt( 1 + 4 * t^2 ) ) / 2;
        z = ( 1 - 1/t ) * x + (1/t) * s;
       
        
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

