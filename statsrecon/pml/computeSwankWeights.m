function w = computeSwankWeights( geom, spectrum , rho, fsoft, fbone, fmetal, sinoMetalMap, metalMaterial  )
% function w = computeSwankWeights( geom, spectrum , rho, fsoft, fbone, fmetal, sinoMetalMap, metalMaterial  )
%
% Meng Wu 
% 2014 at Stanford Univeristy

[ varianceRatios, tableSoft, tableBone, tableMetal,] = computeCompoundPoissonVarianceRatioLookupTable( spectrum , metalMaterial );

%estimate sinogram
ssoft   = forwardProjectMex( fsoft .* rho, geom, 1, 0 );
sbone   = forwardProjectMex( fbone .* rho, geom, 1, 0 );
smetal  = forwardProjectMex( fmetal .* rho, geom, 1, 0 );

w = zeros( size(sinoMetalMap), 'single');
w(sinoMetalMap(:)) = interp3( tableSoft, tableBone, tableMetal, varianceRatios, ...
    ssoft(sinoMetalMap(:)),  sbone(sinoMetalMap(:)),  smetal(sinoMetalMap(:)) );


end