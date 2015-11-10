function steradians = solidAnglePyramidFromDistanceAndBaseLengths(d, a, b)
% steradians = solidAnglePyramidFromDistanceAndBaseLengths(d, a, b)
%    computes the solid angle (in steradians) that is covered by a pyramid
%    with the given distance of center of the pyramid's base rectangle to
%    its apex and the rectangle's base lengths a and b (all in the same
%    length units).
%    See http://en.wikipedia.org/wiki/Solid_angle#Pyramid for the formula.
%
% Copyright (c) 2012 by Andreas Keil, Stanford University.

steradians = 4*atan(a*b/(2*d*sqrt(4*d^2+a^2+b^2)));
