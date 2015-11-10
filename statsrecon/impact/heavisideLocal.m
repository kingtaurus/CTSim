function H = heavisideLocal(phi, epsilon)
% delta = heavisideLocal(phi, epsilon)
% An approximation of the Heaviside function with local support:
% H = 0.5*(1 + phi/epsilon + sin(pi*phi/epsilon)/pi) for |phi| < epsilon,
%     0 for phi <= -epsilon, and
%     1 for phi >= epsilon.
%
% Copyright (c) 2010 by Andreas Keil, Technische Universitaet Muenchen.


H = 0.5*(1 + phi/epsilon + sin(pi*phi/epsilon)/pi);
H(phi <= -epsilon) = 0;
H(phi >= epsilon) = 1;

% if phi <= -epsilon
% 	H = 0;
% elseif phi >= epsilon
% 	H = 1;
% else
% 	H = 0.5*(1 + phi/epsilon + sin(pi*phi/epsilon)/pi);
% end
