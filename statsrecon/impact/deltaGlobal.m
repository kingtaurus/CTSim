function delta = deltaGlobal(phi, epsilon)
% delta = deltaGlobal(phi, epsilon)
% An approximation of the Dirac distribution with global support:
% delta = epsilon/(epsilon^2 + (phi*pi)^2)
%
% Copyright (c) 2010 by Andreas Keil, Technische Universitaet Muenchen.


delta = epsilon./(epsilon^2 + (phi*pi).^2);
