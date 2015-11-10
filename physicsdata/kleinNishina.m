function fKN = kleinNishina(E)
% fKN = kleinNishina(E)
%   returns the value of the Klein-Nishina function for the energy E.
%
% E  photon energy (in keV)
% fKN  value of the Klein-Nishina function (dimensionless)

mec2 = 510.998928; % electron rest energy in keV
alpha = E/mec2;

fKN = (1+alpha)./alpha.^2 .* (2*(1+alpha)./(1+2*alpha) - log(1+2*alpha)./alpha) ...
	+ log(1+2*alpha)./(2*alpha) ...
	- (1+3*alpha)./(1+2*alpha).^2;
