function xf = hilberTransform( x, r, dim )
% function xf = hilberTransform( x, r, dim )
% Hilbert transform with zero padding
% input:
%       x - vector to transform
%       r - ratio of zero padding
%       dim - dimestion of hilber transform if x is a matrx
% output: 
%       xf - hilbert transiform results
%   
% Meng Wu
% 2014

if nargin < 3, dim = 2; end


if isvector(x)
    
    xf = hilber1D( x, r );
    
elseif ismatrix(x)
    
    if dim == 2
        x = x';
    end
    
    xf = zeros( size(x) );
    for i = 1:size(x, 2)
        xf(:,i) = hilber1D( x(:,i), r );
    end
    
    if dim == 2
        xf = xf';
    end
    
else
    error( 'vector or matrix only.');
end


end


function xf = hilber1D( x, r )

n = length(x);
N = 2^ceil( log2( n ) + r );
k = round( ( N - n ) / 2 );

z( k+1: k+n ) = x;
z = imag( hilbert(z) );
xf = z( k+1: k+n );

end