function steradians = solidAnglePyramidFromApexAnlges(a, b)
% steradians = solidAnglePyramidFromApexAnlges(a, b)
%    computes the solid angle (in steradians) that is covered by a pyramid
%    with the given apex angles a and b (in radians).
%    See http://en.wikipedia.org/wiki/Solid_angle#Pyramid for the formula.
%
% Copyright (c) 2012 by Andreas Keil, Stanford University.

steradians = 4*asin(sin(a/2)*sin(b/2));
