function z =  dynamicRangeAdjustment( x, t, eta )
% function z =  dynamicRangeAdjustment( x, t, eta )
%
% Kim, Donghwan, Debashish Pal, Jean-Baptiste Thibault, and Jeffrey A.
% Fessler. “Accelerating Ordered Subsets Image Reconstruction for X-Ray CT
% Using Spatially Nonuniform Optimization Transfer.” IEEE Transactions on
% Medical Imaging 32, no. 11 (November 07, 2013): 1965–78. doi:10.1109/TMI.2013.2266898. 
%
% 2014.4 Meng Wu at Stanford University

z = empiricalCDF( x );
z = max( z.^t, eta );

end

function y = empiricalCDF( x )

a = min( x(:) );
b = mean( x(:) );
c = max( x(:) );

t1 = ( a + b ) / 2;
prob1 = sum( x(:) <= t1 ) / numel(x);

t2 = b;
prob2 = sum( x(:) <= t2 )/ numel(x);

t3 = ( b + c ) / 2;
prob3 = sum( x(:) <= t3 )/ numel(x);

y = zeros( size(x), 'single');
y( x <= t1 ) = prob1 * ( x( x <= t1 ) - a ) / ( t1 - a );
y( x > t1 & x <= t2 ) = ( prob2 - prob1 ) * ( x( x > t1 & x <= t2 ) - t1 ) / ( t2 - t1 ) + prob1;
y( x > t2 & x <= t3 ) = ( prob3 - prob2 ) * ( x( x > t2 & x <= t3 ) - t2 ) / ( t3 - t2 ) + prob2;
y( x > t3 ) = ( 1 - prob3 ) * ( x( x > t3  ) - t3 ) / ( c - t3 ) + prob3;


end


