function w = shortScanModifiedParkerWeight( gammas, beta, Tau, q )
% function m = shortScanParkerWeight( gammas, beta)
%
% Meng Wu at Stanford University
% 2013

if nargin < 4
    q = 0.5;
end

w = zeros( size(gammas) );
for i = 1 : length( gammas )
    
    gamma   = - gammas(i); 
    B       = ( Tau - 2 * gamma );
    b       = q * ( Tau - 2 * gamma );
    bn      = q * ( Tau + 2 * gamma );
    
    w(i)    = 0.5 * ( Sfunction( beta  / b - 0.5 ) ...
                        + Sfunction( ( beta - B ) / b + 0.5 ) ...
                        - Sfunction( ( beta - pi + 2 * gamma ) / bn - 0.5 ) ...
                        - SfunctionS( ( beta - pi - Tau ) / bn + 0.5 ) );
    
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