function w = shortScanParkerWeight( gammas, beta, Tau, nv )
% function m = shortScanParkerWeight( gammas, beta)
%
% Meng Wu at Stanford University
% 2013

if nargin < 4
    nv = 1;
end

w = zeros( size(gammas) );

for i = 1 : length( gammas )
    
    gamma   = - gammas(i);     
    w(i)    = Sfunction( beta / ( Tau - 2 * gamma ) - 0.5 ) ...
        - Sfunction( ( beta - pi + 2 * gamma  ) / ( Tau + 2 * gamma ) - 0.5  ) ;
    
end


if nv > 1
    w = repmat( w, nv, 1 );
end


end


function s = Sfunction( b )

 
if b <= - 0.5
    s = 0;
elseif abs(b) < 0.5
    s = 0.5 * ( 1 + sin( pi *  b) );
else
    s = 1;
end

end