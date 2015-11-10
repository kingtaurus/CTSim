function [ tomoAngles, z ] = tomoAngleCircluarGeometry( geom, scanAngle, z )
% compute tomographic angles for circular geometry tomosynthesis system
%
% Meng Wu at Stanford University
% 2013

if nargin < 3
    z = - 250: 1: 250;
end

% detector size
ldet = geom.detSize(1) * geom.detSpacing(1);

% half scan angle
beta = ( scanAngle ) / 180 * pi / 2;

% source sweep size
lsrc = sin( beta ) * geom.SAD  * 2;

% half fan angle
gamma = atan( ldet / 2 /  geom.SDD );

% compute valid region boundaries
lowerBoundary = - geom.SAD * sin( gamma ) / sin( beta + gamma );
upperBoundary =  geom.SAD * sin( gamma ) / sin( beta - gamma );

if  min(z)  < lowerBoundary || max(z) > upperBoundary
    fprintf('Warning: outside valid region. \n');
end


tomoAngles = zeros( size(z));

for iz = 1:length( z )
    tomoAngles(iz) = asin( lsrc / 2 / (geom.SAD + z(iz) ) ); 
end

tomoAngles = tomoAngles  * 180 / pi * 2;

end