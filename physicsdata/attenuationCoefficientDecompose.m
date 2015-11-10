function [ phi, theta, muEffective, Phi, Theta ] = attenuationCoefficientDecompose( material, spectrum )

energyBinLabels         = spectrum.energyBinLabels;
photonsPerEnergyBin     = spectrum.photonsPerEnergyBin * spectrum.DQE;
energyAverage           = spectrum.energyAverage;


[ muPoly, muEffective] = materialAttenuation( energyBinLabels, material, photonsPerEnergyBin );

Phi = photoElectric(energyBinLabels)./photoElectric(energyAverage);
Theta = kleinNishina(energyBinLabels)./kleinNishina(energyAverage);

b = muPoly;
A = [ Phi, Theta ];

weights = photonsPerEnergyBin;

% apply weights to problem
b = diag(weights) * b;
A = diag(weights) * A;

x = A \ b;
phi = x(1);
theta = x(2);

% compute actual and approximated attenuation at E0
muApprox = Phi * phi + Theta * theta ;

% plot

if false
    figure('Name', ['Attenuation of ', material], 'NumberTitle', 'off');
    zoom on;
    loglog(energyBinLabels, muPoly, 'b-');
    hold on;
    loglog(energyBinLabels, muApprox, 'r-');
    loglog(energyAverage, muEffective, 'bx', 'MarkerSize', 10);
    hold off;
    xlabel('Photon energy (keV)');
    ylabel('Attenuation (1/cm)');
    legend('mu', 'mu_{approx}');
    grid on;
    
end


end