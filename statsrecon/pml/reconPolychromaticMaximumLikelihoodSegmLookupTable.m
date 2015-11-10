function [ u, phis] = reconPolychromaticMaximumLikelihoodSegmLookupTable( sino, ...
    spectrum, geom,  rho, fsoft, fbone, beta, pfun, itnlim,delta, material1, material2 )
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


if nargin < 7
    beta = 1e3;
end

if nargin < 8
    pfun = 'huber';
end

if nargin < 9
    itnlim = 100;
end

if nargin < 10
    delta = 1e-3;
end


if nargin < 11
    material1 = 'Tissue_Soft_ICRU-44';
    material2 = 'Bone_Cortical_ICRU-44';
end


stopCrt         = 1e-5;
showCostFunc    = false;
numos           = 8;

if size( rho, 3) > 10
    k = 2;
else
    k = 1;
end

fprintf('\n\nReconstructed image using polychromatic statistical reconstruction methods: ');
tRecon = tic;

% select projection geometry model
Aos     = @( x, is )forwardProjectMex( x, geom, numos, is, 'proj,dd' );
Atos    = @( x, is )backProjectMex( x, geom, numos, is, 'back,dd' );

A       = @( x )forwardProjectMex( x, geom, 1, 0, 'proj,dd' );
At      = @( x )backProjectMex( x, geom, 1, 0, 'back,dd' );
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

% tabulate Ybar and nablaY
[ Ybar, nablaYsoft, nablaYbone, tableSoft, tableBone, muSoft, muBone] = computePolychromaticAttLookupTable( spectrum, material1, material2 );

% precompute curvature

E = muSoft * fsoft + muBone * fbone;
%curv = ( muSoft^2+ muBone^2 * fbone ) .* At( A( ones(size(rho), 'single' ) ) .* sino );
curv = E .*  At( A( E  ) .* sino );
u = E .* rho ;

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
        ssoft = Aos( fsoft .* rho, isub );
        sbone = Aos( fbone .* rho, isub );
        
        ssoft(ssoft > 60 ) = 60;
        sbone(sbone > 20 ) = 20;
        
        y   = interp2( tableSoft, tableBone, Ybar, ssoft, sbone);
        dysoft = interp2( tableSoft, tableBone, nablaYsoft, ssoft, sbone);
        dybone = interp2( tableSoft, tableBone, nablaYbone, ssoft, sbone);
        
        % intensity correction if use bowtie filter
        if spectrum.useBowtie
            
            if ndims( y ) == 2
                
                for iv = 1 : size(y, 2)
                    y(:, iv) = y(:, iv) .* spectrum.flatFieldRatio;
                    dysoft(:, iv) = dysoft(:, iv) .* spectrum.flatFieldRatio;
                    dybone(:, iv) = dybone(:, iv) .* spectrum.flatFieldRatio;
                end
                
                
            else
                
                for iv = 1 : size(y, 3)
                    y(:, :, iv) = y(:, :, iv) .* spectrum.flatFieldRatio;
                    dysoft(:, :, iv) = dysoft(:, :, iv) .* spectrum.flatFieldRatio;
                    dybone(:, :, iv) = dybone(:, :, iv) .* spectrum.flatFieldRatio;
                end
                
            end
        end
        
        
        %compute gradient
        p = ( Os(sino, isub) ./ y ) - 1;
        N =  fsoft .* Atos(  p.* dysoft, isub )  +  fbone .* Atos(  p.* dybone, isub );
        
        
        %piecewise smooth regulator
        numerator =  - N + beta * E .* S(u);
        denominator = curv + beta * ( E ) .*  T(u);
        
        %update the sinogram
        rho = rho - numerator ./ denominator;
        
        rho = extendVoi( rho, k );
        rho( rho < 0 ) = 0;
        rho( isnan(rho) ) = 0;
        rho( isinf(rho) ) = 0;
        
        if isub == 1
            phi = 0;
        end
        yy = Os(sino, isub);
        phi = phi + sum( ( yy(:) - y(:) ).^2 ./ yy(:) );
        
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

