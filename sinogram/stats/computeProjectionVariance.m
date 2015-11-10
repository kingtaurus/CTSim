function sinoVariance = computeProjectionVariance( sinoPC, elecNoiseVar, gain )
% sinoVariance = computeProjectionVariance( sinoPC, electronicNoiseStd )
% Compute the variance of the logged sinogram
%   
% Meng Wu @ Stanford University
% 2014

if nargin < 3
    gain = 1;
end

a = floor( log10( min(sinoPC(:) ) ));
b = ceil( log10( max(sinoPC(:) ) ));

[photonCounts , variance ] = generateLogProjectionVarianceLookupTable( sqrt( elecNoiseVar) , a, b );

% figure;
% loglog(photonCounts, variance);


sinoVariance = interp1( photonCounts,  variance,  sinoPC) / gain;