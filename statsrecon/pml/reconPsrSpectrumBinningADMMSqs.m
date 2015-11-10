function [ x, phis] = reconPsrSpectrumBinningADMMSqs( sino, ...
    BinSizes, Phis, Thetas, geom,  rho, fphi, ftheta, phisFiltraion, thetaFiltrations, ...
    beta, pfun, itnlim, delta, numos )
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



if nargin < 12
    pfun = 'huber';
end

if nargin < 13
    itnlim = 100;
end

if nargin < 14
    delta = 1e-3;
end

if nargin < 15
    numos   = 8;
end

stopCrt         = 1e-5;
showCostFunc    = false;
numBins         = length(BinSizes);

fprintf('\n\n Polychromatic Statistical Reconstruction using %d - Spectrum Binning using ADMM SQS: ', numBins);
tRecon = tic;

% select projection geometry model
Aos     = @( x, is )forwardProjectMex( x, geom, numos, is );
Atos    = @( x, is )backProjectMex( x, geom, numos, is );
A2os     = @( x1, x2, is )forwardProject2Mex( x1, x2, geom, numos, is, 'proj,dd' );
A2tos    = @( x1, x2, is )backProject2Mex( x1, x2, geom, numos, is, 'back,dd' );
A       = @( x )forwardProjectMex( x, geom, 1, 0 );
At      = @( x )backProjectMex( x, geom, 1, 0 );
Os      = @( x, is)orderedSubset(x, numos , is);

C   = @( x )piecewiseDifference( x, geom.reconSize, 1 );
Ct  = @( x )piecewiseDifference( x, geom.reconSize, 2 );
m = piecewiseDifference( rho, geom.reconSize, 3 );

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
E =  fphi +  ftheta;
x = E .* rho ;
curv = E.^2 .* At( A( ones( size( rho ), 'single' ) ) .* sino );

fprintf('\nbeta  = %11.2e', beta );
fprintf('\ndelta = %10.3f', delta );
fprintf('\titnlim = %10g', itnlim);
hdg1 = '   itn       x(0)          PHI(x)        b*R(x)        RMSD';
fprintf('\n%s'      , hdg1  );


phis = zeros(itnlim, 1);
rmsds = zeros(itnlim, 1);

gamma = 1;
u = C(E .* rho) ;
d =  zeros( size(u), 'single' );
% start iterations
for itn = 1 : itnlim
    
    xold = x;
    
    % step 1: rho^{(k+1)} = g( rho ) + gamma / 2 * || C * E * rho - u^{(k+1)} - d^{(k+1)} ||^2
    for isub = 1:numos
        
        %estimate sinogram
        [sphi, stheta] = A2os( fphi .* rho, ftheta .* rho, isub );
        
        % intensity correction if use bowtie filter
        if min( phisFiltraion(:) ) > 0
            if ndims( sphi ) == 2
                for iv = 1 : size(sphi, 2)
                    sphi(:, iv) = sphi(:, iv) + phisFiltraion;
                    stheta(:, iv) = stheta(:, iv) + thetaFiltrations;
                end
            else
                for iv = 1 : size(sphi, 3)
                    sphi(:, :, iv) = sphi(:, :, iv) + phisFiltraion;
                    stheta(:, :, iv) = stheta(:, :, iv) + thetaFiltrations;
                end
            end
        end
        
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
            
            temp    = BinSizes(ib) * exp( - sphi * Phis(ib) - stheta * Thetas(ib) );
            y       = y + temp;
            dyPhi  = dyPhi - Phis(ib) * temp;
            dyTheta  = dyTheta - Thetas(ib) * temp;
        end
        
        
        %compute gradient
        p = ( Os(sino, isub) ./ y ) - 1;
        
        [lphi, ltheta ] = A2tos( p.* dyPhi,  p.* dyTheta, isub );
        N =  fphi .* lphi +  ftheta .* ltheta;
        
        %update the sinogram
        b = beta * max( 2 - 2 * itn / itnlim , 1 );
        rho = rho - (- N + gamma * b  * E.* Ct( C( E .* rho ) - u - d  )) ./ (curv + gamma * b * E.* Ct( C( E ) ));
        
        rho( rho < 0 ) = 0;
        rho( isnan(rho) ) = 0;
        rho( isinf(rho) ) = 0;
        
        if isub == 1
            phi = 0;
        end
        yy = Os(sino, isub);
        phi = phi + sum( ( yy(:) - y(:) ).^2 ./ yy(:) );
        
        x = E .* rho;
        s = C( E .* rho ) - d;
        %step 2: u^{(k+1)} = argmin h(u) + gamma / 2 || C * E * rho^{(k+1)} - u - d^{(k)} ||^2
        for t = 1:10
            u = u - ( m .* dHuber( u, delta ) + gamma * ( u - s ) ) ./ ( m .* cHuber( u, delta ) + gamma );
        end
        
        % step 3: d^{(k+1)} = d^{(k)} - C * E * rho^{(k)} + u^{(k+1)}
        d = u - s;
        
    end
    
    
    %print out iteration parameters
    phi = phi +  beta * R(x);
    rmsd = sqrt( mean( (x(:) - xold(:)).^2 ) );
    
    phis(itn) = phi;
    rmsds(itn) = rmsd;
    
    
    prnt = 0;
    if itn   <= 5       , prnt = 1; end
    if itn   >= itnlim-5, prnt = 1; end
    if rem(itn,5) == 0  , prnt = 1; end
    
    if prnt
        fprintf('\n%6g %13.3e'  , itn  , x(round(end/2),round(end/2),ceil(end/2)));
        fprintf(' %13.3e %13.3e %13.3e', phi  , beta * R(x), rmsd);
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
    ylabel '\phi(x)'; xlabel 'No. of iteration';
    
end

end


function d = dHuber( d, delta)
d( d > delta )  = delta;
d( d < -delta ) = - delta;
end

function d = cHuber( d, delta)
qmap = ( abs(d) <= delta );
d( qmap ) = 1;
d( ~qmap ) = 0;
end

