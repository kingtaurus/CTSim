function [photonCounts , variance ] = generateLogProjectionVarianceLookupTable( electronicNoiseStd, a, b )

if nargin < 2
    a = 2;
    b = 8;
end

N = 100000;
M = 200;
photonCounts = logspace(a, b, M);

mu = zeros(1,M);
variance = zeros(1,M);


for i = 1:M
    
    if photonCounts(i) < max( electronicNoiseStd^2 * 100, 10000)
        x = poissrnd(photonCounts(i), [N 1]) + electronicNoiseStd * randn([N, 1]);
        mu(i) = mean( log(x));
        variance(i) = var( log(x));
    else
        mu(i) = log( photonCounts(i) );
        variance(i) = 1 / photonCounts(i);
    end
end