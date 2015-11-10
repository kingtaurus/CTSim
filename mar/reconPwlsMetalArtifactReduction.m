function [img, phis] = reconPwlsMetalArtifactReduction( sino, spectrum, geom, weights, sinoMapUnknown,  beta, delta, itnlim )
% Penalized-Weight-Least-Squares recosntruction for single energy sinograms.
% input:
%       sino    - photon counts sinogram
%       spectrum - X-ray spectrum infomation
%       geom    - system geometry
%       beta    - roughness penalty weight (default 0.01)
%       delta   - parameters for huber function (default 0.1)
%       cstop   - cost function stopping criteria
%       nitn    - maximum number of iterations
%
% output:
%       img     - reconstruction result
%       phis    - cost function values
%
% Based "Statistical image reconstruction for polyenergetic X-ray computed
% tomography" by Elbakri, Idris A and Fessler, Jeffrey A
%
% 2012.10 Meng Wu at Stanford University

if nargin < 6
    beta = 1e5;
end

if nargin < 7
    delta = 0.1;
end

if nargin < 8
    itnlim = 50;
end

stopCrt         = 1e-4;
showCostFunc    = false;
M               = 16;

% select projection geometry model
Aos     = @( x, is )forwardProjectMex( x, geom, M, is );
Atos    = @( x, is )backProjectMex( x, geom, M, is );
Os      = @( x, is )orderedSubset( x, M , is);
A       = @( x )forwardProjectMex( x, geom );
At      = @( x )backProjectMex( x, geom );
R       = @( x )huberPenalty( x, 0, delta );
S       = @( x )huberPenalty( x, 1, delta );
T       = @( x )huberPenalty( x, 2, delta );
% R       = @( x )quadPenalty( x, 0 );
% S       = @( x )quadPenalty( x, 1 );
% T       = @( x )quadPenalty( x, 2 );

% estimate line intergrals
[u, l ] = reconInterpFBP( sino, geom, spectrum, sinoMapUnknown  );

w = weights;
w(sinoMapUnknown)  = mean( weights(:) / 10 );
a = A( ones( size(u), 'single' ) );
precom = At( w .* a );
w(sinoMapUnknown) = 0;

fprintf('\nbeta  = %11.2e', beta );
fprintf('\ndelta = %10.3f', delta );
fprintf('\titnlim = %10g', itnlim);
hdg1 = '   itn       x(0)          PHI(u)        b*R(u)';

fprintf('\n%s'      , hdg1  );

phis = zeros(itnlim, 1);

for itn = 1 : itnlim
    
    for isub = 1 : M
        
        d = Aos( u, isub )- Os(l, isub);
        
        numerator = Atos( Os(w, isub).* d, isub ) + beta * S(u);
        
        denumerator = precom + beta * T(u);
        
        u = u - numerator ./ denumerator;
        
        u( u < 0 ) = 0;
        u( isnan(u) ) = 0;
        u( isinf(u) ) = 0;
        
        
        if isub == 1
            phi = 0;
        end
        
        temp = Os(w, isub ).* ( d.^2 );
         phi = phi + sum( temp(:) ) / 2;
    
    end
    
    phi = phi + beta * R(u);
    phis(itn) = phi;
    
    prnt = 0;
    if itn   <= 5       , prnt = 1; end
    if itn   >= itnlim-5, prnt = 1; end
    if rem(itn,10) == 0  , prnt = 1; end
    
    if prnt
        fprintf('\n%6g %13.3e'  , itn  , u(round(end/2),round(end/2),ceil(end/2)));
        fprintf(' %13.3e %13.3e', phi  , beta * R(u))
    end
    
    if itn > 1 && abs( phi - phis(itn-1) ) / phi < stopCrt
        break;
    end
       
    
end

img = u;

fprintf('\n\n');

if showCostFunc
    figure('Name', 'cost function values', 'NumberTitle', 'off');
    semilogy( 1:itn, phis(1:itn) );
    ylabel '\phi(u)'; xlabel 'No. of iteration';
    
end


end


