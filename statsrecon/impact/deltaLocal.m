function delta = deltaLocal(phi, epsilon)
% delta = deltaLocal(phi, epsilon)
% An approximation of the Dirac distribution with local support:
% delta = (1 + cos(pi*phi/epsilon))/(2*epsilon) for |phi| < epsilon and
%         0 otherwise.
%
% Copyright (c) 2010 by Andreas Keil, Technische Universitaet Muenchen.


delta = (1 + cos(pi*phi/epsilon))/(2*epsilon);
delta(abs(phi) >= epsilon) = 0;

% if abs(phi) < epsilon
% 	delta = (1 + cos(pi*phi/epsilon))/(2*epsilon);
% else
% 	delta = 0;
% end
