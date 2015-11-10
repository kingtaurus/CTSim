function xf = filterFreq2D( x, H )
% function filterFreq( x, H  )
% filtering signal with H filter
% input:
%       x - vector to transform
%       H - filter in frequency domain
% output:
%       y - hilbert transiform results
%
% Meng Wu
% 2014
[nv, nu] = size(x);
[mv, mu] = size(H);

y = zeros( [mv, mu] );
% filter projection slice
y(1:nv,1:nu) = x;
pFFT = (fft2(y)); % FFT of projection
pFFT = pFFT .* H; % filter in frequency domain
y = real(ifft2(pFFT)); % use real part of inverse FFT
xf = y(1:nv,1:nu); % truncate


end