function intensity = simulatePhotonCountingData( phan, geom, spectrum, sinosDir, usePolyEnergies, addCountingNoise )

if nargin < 5
    usePolyEnergies = true;
end


if nargin < 6
    addCountingNoise = true;
end


fprintf('The name has been changed to simulateSimplePoissonData.m. \n');

intensity = simulateSimplePoissonData( phan, geom, spectrum, sinosDir, usePolyEnergies, addCountingNoise );




end