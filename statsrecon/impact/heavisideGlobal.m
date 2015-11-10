function H = heavisideGlobal(phi, epsilon)
% delta = heavisideGlobal(phi, epsilon)
% An approximation of the Heaviside function with global support:
% H = 0.5 + atan(pi*phi/epsilon)/pi
%
% Copyright (c) 2010 by Andreas Keil, Technische Universitaet Muenchen.


H = 0.5 + atan(pi*phi/epsilon)/pi;
