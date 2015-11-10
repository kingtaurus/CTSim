function yi = interp1geom(x, y, xi, varargin)
% yi = interp1geom(x, y, xi, varargin)
%   performs interpolations just like the built-in interp1 function with
%   the difference that all values are interpolated geometrically instead
%   of arithmetically, i.e., interp1geom(...) = exp(interp1(log(...))).
%
% Copyright (c) 2012 by Andreas Keil, Stanford University.

yi = exp(interp1(log(x), log(y), log(xi), varargin{:}));
