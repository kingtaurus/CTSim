function [phan, map ] = loadMaterialsDensityPhantom( p )
% Load materials and density defined Phantom
%
% Meng Wu at Stanford University
% 2014 - 2015

% read phantom material data file
[phan, map] = loadMaterialsPhantom( p );

% read phantom density data file
fprintf( '\t and load phantom density map ... ' );

phantomDensities = readMetaImage( [ p.Phantom.materialsFileName '-density.mhd' ] );

% 2D case
if length(phan.reconSize) == 2
    phantomDensities = phantomDensities(:,:,ceil(end/2));
end

phan.phantomDensities       = phantomDensities;
phan.useDensity             = 1;


fprintf('done.\n');

end