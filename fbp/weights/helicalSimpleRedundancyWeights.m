function weights = helicalSimpleRedundancyWeights( gammas, beta, Delta )

a = Delta / 2;
b = Delta / 2 - 0.5;

lambda = beta - a;
lambdas = lambda + (-3:3) * pi;
lambda1 = abs( lambdas );

weights = zeros( size(gammas) );

for i = 1:length(gammas)
    
    lambda2 = abs(  lambdas + [-2 0 -2 0 -2 0 -2] * gammas(i) );

    c = ( sin( pi * ( a - lambda2 ) ) ) .^2;
    c( lambda2 <= b ) = 1;
    c( lambda2 > a ) = 0;

    
    weights(i) =  c(4) / sum(c);
    
end

weights = conv( weights, [0.25 0.5 0.25], 'same' );

end

