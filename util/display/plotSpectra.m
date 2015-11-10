%function plotSpectra( spectrum )

figure;

[energyBinLabels1, ~, photonsPerEnergyBin1, ~, ~, ~] = ...
    readSpectrum( 'spectrum_2.5MeVp_1mm_Tungsten.txt', 'per bin', 1);


[energyBinLabels2, ~, photonsPerEnergyBin2, ~, ~, ~] = ...
    readSpectrum( 'spectrum_6MeVp_WCuTarget_ClinacFilter.txt', 'per bin', 1);

plot( energyBinLabels2 / 1000, photonsPerEnergyBin2 / max( photonsPerEnergyBin2),  '--', 'lineWidth', 2, 'Color', 'k'  );
hold on;
plot( energyBinLabels1 / 1000, photonsPerEnergyBin1 / max( photonsPerEnergyBin1), 'lineWidth', 2, 'Color', 'k'  ); 

xlabel( 'MeV', 'fontSize', 20);
legend( '6 MVp', '2.5 MVp'  );

set(gca,'FontSize',20);
axis tight; grid on;
 saveas(gcf, 'spec.eps');



% plot( spectrum.energyBinLabels / 1000, spectrum.photonsPerEnergyBin / max( spectrum.photonsPerEnergyBin), 'lineWidth', 2, 'Color', 'k'  );
% xlabel( 'MeV', 'fontSize', 20);
% set(gca,'FontSize',20);
% axis tight;
%  saveas(gcf, 'spec.eps');