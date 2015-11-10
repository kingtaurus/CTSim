function [ u, phis] = reconPolychromaticMaximumLikelihoodSpectrumBinning( sino, ...
    binPhotonCounts, Phis, Thetas, geom,  rho, fphi, ftheta, beta, pfun, itnlim,delta )
% Statistical reconstruction for polyenergetic CT for single energy sinograms.
% input:
%       sino    - photon counts sinogram
%       spectrum - X-ray spectrum infomation
%       geom    - system geometry
%       beta    - roughness penalty weight (default 0.01)
%       delta   - parameters for huber function (default 0.1)
%       itnlim    - maximum number of iterations
%
% output:
%       img     - reconstruction result
%       phis    - cost function values
%
% Based "Statistical image reconstruction for polyenergetic X-ray computed
% tomography" by Elbakri, Idris A and Fessler, Jeffrey A
%
% 2012.10 Meng Wu at Stanford University


if nargin < 9
    beta = 1e3;
end

if nargin < 10
    pfun = 'huber';
end

if nargin < 11
    itnlim = 100;
end

if nargin < 12
    delta = 1e-3;
end

stopCrt         = 1e-5;
showCostFunc    = false;
numos           = 36;

if size( rho, 3) > 10
    k = 2;
else
    k = 1;
end

fprintf('\n\n Polychromatic Maximum Likelihood Reconstruction using %d - Spectrum Binning: ', length(Phis));
tRecon = tic;

% select projection geometry model
Aos     = @( x, is )forwardProjectMex( x, geom, numos, is );
Atos    = @( x, is )backProjectMex( x, geom, numos, is );
A2os     = @( x1, x2, is )forwardProject2Mex( x1, x2, geom, numos, is, 'proj,dd' );
A2tos    = @( x1, x2, is )backProject2Mex( x1, x2, geom, numos, is, 'back,dd' );
A       = @( x )forwardProjectMex( x, geom, 1, 0 );
At      = @( x )backProjectMex( x, geom, 1, 0 );
Os      = @( x, is)orderedSubset(x, numos , is);

% select penalty function
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

% precompute curvature
curv = ( fphi + ftheta ) .* At( A( ones(size(rho), 'single' ) ) .* sino );
E =  fphi +  ftheta;
u = E .* rho ;

numBins = length(binPhotonCounts);

fprintf('\nbeta  = %11.2e', beta );
fprintf('\ndelta = %10.3f', delta );
fprintf('\titnlim = %10g', itnlim);
hdg1 = '   itn      rho(0)          PHI(u)        b*R(u)';
fprintf('\n%s'      , hdg1  );



phis = zeros(itnlim, 1);
rmsds = zeros(itnlim, 1);

% start iterations
for itn = 1 : itnlim
    
    uold = u;
    
    for isub = 1:numos
        
        %estimate sinogram
        [sphi, stheta] = A2os( fphi .* rho, ftheta .* rho, isub );
        
        
        sphi(sphi > 60 ) = 60;
        stheta(stheta > 30 ) = 30;
        
        if itn > 1 &&   mod(geom.noViews, numos) == 0
            y(:) = 0;
            dyPhi(:) = 0;
            dyTheta(:) = 0;
            
        else
            y       = zeros( size(sphi), 'single' );
            dyPhi  = zeros( size(sphi), 'single' );
            dyTheta  = zeros( size(sphi), 'single' );
            
        end
        
        
        for ib = 1: numBins
            
            temp    = binPhotonCounts(ib) * exp( - sphi * Phis(ib) - stheta * Thetas(ib) );
            y       = y + temp;
            dyPhi  = dyPhi - Phis(ib) * temp;
            dyTheta  = dyTheta - Thetas(ib) * temp;
        end
        
        
        %compute gradient
        p = ( Os(sino, isub) ./ y ) - 1;
        
        [lphi, ltheta ] = A2tos( p.* dyPhi,  p.* dyTheta, isub );
        N =  fphi .* lphi +  ftheta .* ltheta;
        
        %piecewise smooth regulator
        numerator =  - N + beta .* S(u);
        denominator = curv + beta * E .*  T(u);
        
        %update the sinogram
        rho = rho - 1.5 * numerator ./ denominator;
        
        rho = extendVoi( rho, k );
        rho( rho < 0 ) = 0;
        rho( isnan(rho) ) = 0;
        rho( isinf(rho) ) = 0;
        
        if isub == 1
            phi = 0;
        end
        phi = phi + sum( sum( sum( (log( p + 1) ).^2 .* Os(sino, isub) )) ) / 2;
        
    end
    
    u = E .* rho;
    %print out iteration parameters
    phi = phi +  beta * R(u);
    rmsd = sqrt( mean( (u(:) - uold(:)).^2 ) );
    
    phis(itn) = phi;
    rmsds(itn) = rmsd;
    
    
    prnt = 0;
    if itn   <= 5       , prnt = 1; end
    if itn   >= itnlim-5, prnt = 1; end
    if rem(itn,5) == 0  , prnt = 1; end
    
    if prnt
        fprintf('\n%6g %13.3e'  , itn  , u(round(end/2),round(end/2),ceil(end/2)));
        fprintf(' %13.3e %13.3e %13.3e', phi  , beta * R(u), rmsd);
    end
    
    if itn > 5 && rmsd < stopCrt
        break;
    end
    
    
end

tRecon = toc(tRecon);
fprintf('\nDone in %dmin %0ds.\n', floor(tRecon/60), round(mod(tRecon, 60)));

if showCostFunc
    figure('Name', 'cost function values', 'NumberTitle', 'off');
    semilogy( 1:itn, phis(1:itn) );
    ylabel '\phi(u)'; xlabel 'No. of iteration';
    
end

end

