function [ u, phis] = reconIterativeSinogramPatch(  sinoPCKeV, sinoPCMeV, ...
    spectrumKeV, spectrumMeV, geomKeV, geomMeV,  beta, itnlim, delta )
% PWLS reconstruction for polyenergetic CT for dual energy sinograms.
% input:
% output:
%       img     - reconstruction result
%       phis    - cost function values
%
% 2012.10 Meng Wu at Stanford University
%
% This a non convex version of pwls. The problem come from the fact that
% the conversion between the keV and MeV image are not convex or concave.
% Moreover, the lack of the information of metal in keV image create additional problem.


if nargin < 7
    beta = 1e3;
end

if nargin < 8
    pfun = 'huber';
end

if nargin < 9
    itnlim = 50;
end

stopCrt         = 1e-4;
showCostFunc    = false;
M               = 2;
delta           = 1e-4;

% select projection geometry model
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

fprintf('Reconstructing with keV/MeV data using iterative sinogram patching... \n');
tRecon = tic;

% select projection geometry model
rk = spectrumKeV.backgroundEvents;
rm = spectrumMeV.backgroundEvents;
sinoKeV = computeSinogramAttunation( sinoPCKeV , spectrumKeV );
sinoMeV = computeSinogramAttunation( sinoPCMeV , spectrumMeV );


[sinoMeV, sinoMetalMap, sinoOverlapMap] = collimateMevSinogram(  sinoKeV, sinoMeV, ...
    geomKeV, geomMeV, spectrumKeV, 300 );

% compute initial image use analytical method
[u, l ]= reconFBPSinogramDEPatch1( sinoKeV, sinoMeV, geomKeV, ...
    sinoMetalMap, sinoOverlapMap, 0, spectrumKeV, spectrumMeV );

sinoPCMeV = convertSinogramGeometry(sinoPCMeV, geomMeV, geomKeV );

% compute the weights for two energies
w = ( sinoPCKeV - rk ).^2 ./ sinoPCKeV;
w( sinoPCKeV <= rk) = 0;
w( sinoMetalMap ) = ( sinoPCMeV( sinoMetalMap ) - rm ).^2 ./ sinoPCMeV( sinoMetalMap );

clear sinoPCKeV sinoPCMeV;

% precompute curvature
a = A( ones( size(u), 'single'));
precom = At(  w.* a );

fprintf('\nbeta  = %11.2e', beta );
fprintf('\ndelta = %10.3f', delta );
fprintf('\titnlim = %10g', itnlim);
hdg1 = '   itn      pho(0)          PHI(u)        b*R(u)';

fprintf('\n%s'      , hdg1  );

phis = zeros(itnlim, 1);
for itn = 1 : itnlim
    
    for isub = 1 : M
        
        d = Aos( u, isub ) - Os(l, isub);
        
        numerator = Atos( Os(w, isub).* d, isub ) + beta * S(u) ;
        denumerator = precom + beta * T(u);
        
        u = u -  numerator ./ denumerator;
        
        u( u < 0 ) = 0;
        u( isnan(u) ) = 0;
        u( isinf(u) ) = 0;
        
        if isub == 1
            phi = 0;
        end
        
        temp =  Os(w, isub).* d.^2 ;
        phi = phi + sum( temp(:) ) / 2;
        
    end
    
    phis(itn) = phi +  beta * R(u);
    
    prnt = 0;
    if itn   <= 5           , prnt = 1; end
    if itn   >= itnlim-5    , prnt = 1; end
    if rem(itn,10) == 0     , prnt = 1; end
    
    if prnt
        fprintf('\n%6g %13.3e'  , itn  , u(round(end/2),round(end/2),ceil(end/2)));
        fprintf(' %13.3e %13.3e ', phi  ,  beta * R(u));
    end
    
    if itn > itnlim/2 && abs( phi - phis(itn-1) ) / phi < stopCrt
        break;
    end
    
    % repatch the sinogram
    if mod(itn, 3) == 0
        
        l = 0.5 * ( l + patchSinogramDualEnergies2( u, sinoMeV, sinoKeV, ...
            sinoMetalMap, sinoOverlapMap, geomKeV, spectrumKeV, spectrumMeV ));
    end
    
end


tRecon = toc(tRecon);
fprintf('Done in %dmin %0ds.\n', floor(tRecon/60), round(mod(tRecon, 60)));

if showCostFunc
    figure('Name', 'cost function values', 'NumberTitle', 'off');
    semilogy( 1:itn, phis(1:itn) );
    ylabel '\phi(u)'; xlabel 'No. of iteration';
    
end


end

