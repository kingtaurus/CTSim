function [A, At, Aos, Atos, Os ] = loadPojectors( geom, numos )


if nargin < 2
    numos = 1;
end

Aos     = @( x, is )forwardProjectMex( x, geom, numos, is, 'proj,tf' );
Atos    = @( x, is )backProjectMex( x, geom, numos, is, 'back,tf' );
A       = @( x )forwardProjectMex( x, geom, 1, 0, 'proj,tf' );
At      = @( x )backProjectMex( x, geom, 1, 0, 'back,tf' );
Os      = @( x, is)orderedSubset(x, numos , is);

end