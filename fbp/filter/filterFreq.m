function xf = filterFreq( x, H, dim )
% function filterFreq( x, H, dim  )
% filtering signal with H filter
% input:
%       x - vector or matrix to transform
%       H - filter in frequency domain
%       dim - dimestion of filtering if x is a matrx
% output:
%       xf - filtered signal
%
% Meng Wu
% 2014

if nargin < 3, dim = 2; end
H = H(:);

if isvector(x)
    
    xf = filt1D( x(:), H);
    
elseif ismatrix(x)
    
    if dim == 2
        x = x';
    end
    
    xf = zeros( size(x) );
    for i = 1:size(x, 2)
        xf(:,i) = filt1D( x(:,i), H );
    end
    
    if dim == 2
        xf = xf';
    end
    
else
    error( 'vector or matrix only.');
end

end

function yf = filt1D( y, H )

n = length(y);
% filter projection slice
y(length(H)) = 0; % pad to filter size with zeros
pFFT = (fft(y)); % FFT of projection
pFFT = pFFT .* H; % filter in frequency domain
yf = real(ifft(pFFT)); % use real part of inverse FFT
yf(n+1:end) = []; % truncate

end