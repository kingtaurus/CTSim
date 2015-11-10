function photons = hDetectorEnergyResponseFunction(photons, energy)

scintillatorThickness = 1.5; % cm
scintillatorMaterial = 'CWO'; % 
detectorGeometryLoss = 0.90; %

mu = materialAttenuation( energy, scintillatorMaterial );

photons = detectorGeometryLoss * photons .* ( 1 - exp( - mu .* scintillatorThickness )) ;


end