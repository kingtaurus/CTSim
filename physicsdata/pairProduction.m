function fPP = pairProduction(E)
% fPP = pairProduction(E)
%   returns the value of the energy-dependence of the pair production
%   effect (which approximated by ln(E/E_rest)).
%
% E  photon energy (in keV)
% fPP  value of the Klein-Nishina function (dimensionless)

E_rest = 1.022e3;

fPP = (E > E_rest) .* log(E./E_rest);
