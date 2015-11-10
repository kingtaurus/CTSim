function fPE = photoElectric(E)
% fPE = photoElectric(E)
%   returns the value of Andreas Keil's approximation to the
%   energy-dependence of the photoelectric effect (which is often said to
%   be E^(-3) but more closely approximated by E^(-2.6)).
%
% E  photon energy (in keV)
% fPE  value of the Klein-Nishina function (dimensionless)

fPE = E.^(-3);
