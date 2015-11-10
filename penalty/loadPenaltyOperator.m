function [R, S, T, H ]  = loadPenaltyOperator( pfunc, delta )
% function [R, S, T, H ]  = loadPenaltyOperator( pfunc, delta )
% load operators for penalty function
% inputs:
%       pfunc - name of penalty function
%               'huber'  huber function
%               'quad'   quadratic function
%               'hyper'  hyperbolic function
%               'aniso'  anisotropic quadratic function with gaussian weighting
%               'isohuber' isotropic huber function
%               'isotv'    isotropic total variation
%       delta - one paramater (may not be used )
% outputs:
%       R   - penalty function score
%       S   - first derivative
%       T   - second derivative
%       H   - optimal huber curvature
%
% Meng Wu at Stanford University
% 2014

if strcmpi(pfunc, 'huber')
    R       = @( x )huberPenalty( x, 0, delta );
    S       = @( x )huberPenalty( x, 1, delta );
    T       = @( x )huberPenalty( x, 2, delta );
elseif strcmpi(pfunc, 'quad')
    R       = @( x )quadPenalty( x, 0 );
    S       = @( x )quadPenalty( x, 1 );
    T       = @( x )quadPenalty( x, 2 );
elseif strcmpi(pfunc, 'hyper')
    R       = @( x )hyperbolaPenalty( x, 0, delta );
    S       = @( x )hyperbolaPenalty( x, 1, delta );
    T       = @( x )hyperbolaPenalty( x, 2, delta );
elseif strcmpi(pfunc, 'aniso')
    R       = @( x )anisotropicPenalty( x, 0, delta );
    S       = @( x )anisotropicPenalty( x, 1, delta );
    T       = @( x )anisotropicPenalty( x, 2, delta );
elseif strcmpi(pfunc, 'isohuber')
    R       = @( x )huberIsoPenalty( x, 0, delta );
    S       = @( x )huberIsoPenalty( x, 1, delta );
    T       = @( x )huberIsoPenalty( x, 2, delta );
elseif strcmpi(pfunc, 'isotv')
    R       = @( x )totalVariationIsoPenalty( x, 0 );
    S       = @( x )totalVariationIsoPenalty( x, 1 );
    T       = @( x )totalVariationIsoPenalty( x, 2 );
elseif strcmpi(pfunc, 'tv')
    R       = @( x )totalVariationPenalty( x, 0 );
    S       = @( x )totalVariationPenalty( x, 1 );
    T       = @( x )totalVariationPenalty( x, 2 );
else
    error('unknow penalty function');
end
H       = @( x )huberPenalty( x, 3, delta );


end