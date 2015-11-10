function [ doseTotal, doseKeV, doseMeV] = computeDose( p, spectrumKeV, spectrumMeV, outputDir )
% Compute Resulting Dose
%
% Copyright (c) 2010-2012 by Andreas Keil, Stanford University.
% Modified by Meng Wu at 2012.9


% determine directory with dose simulation results
doseDir = [ p.Phantom.materialsFileName(1:end-4) '/'];

% simulation parameters
photonsPerSteradianSimKeV = 1e8*180/pi; % #photons/sr used for the dose simulation
doseSimKeV = 1e-4; % Gy as computed in the dose simulation (dose values sampled and averaged from several points for the simulated #photons/sr)
noViewsSimKeV = 650; % #views used for the dose simulation

% simulation parameters
photonsPerSteradianSimMeV = 1e8*180/pi; % #photons/sr used for the dose simulation
doseSimMeV = 0.00042; % Gy as computed in the dose simulation (dose values sampled adjacent to the fillings for the simulated #photons/sr)
noViewsSimMeV = 650; % #views used for the dose simulation


% load keV dose image
doseKeVSim = readMetaImage([doseDir p.Spectra.spectrumKeV '.mhd']);
doseKeVSim = imresize(doseKeVSim, p.Reconstruction.size, 'box', 'Antialiasing', true);
doseScalingKeV = spectrumKeV.photonsPerSteradian/photonsPerSteradianSimKeV * p.Geometries.keV.noViews/noViewsSimKeV; % = p.Spectra.doseLimitTotalKeV/doseSimKeV
doseKeV = doseKeVSim * doseScalingKeV;

% load MeV dose image
doseMeVSim = readMetaImage([doseDir p.Spectra.spectrumMeV '.mhd']);
doseMeVSim = imresize(doseMeVSim, p.Reconstruction.size, 'box', 'Antialiasing', true);
doseScalingMeV = spectrumMeV.photonsPerSteradian/photonsPerSteradianSimMeV * p.Geometries.MeV.noViews/noViewsSimMeV; % = p.Spectra.doseLimitTotalMeV/doseSimMeV
doseMeV = doseMeVSim * doseScalingMeV;

% compute keV+MeV dose
doseTotal = doseKeV + doseMeV;


if nargin == 4
    figure('Name', 'Total dose', 'NumberTitle', 'off');
    zoom on;
    %	imshow(doseTotal', stretchlim(doseTotal, [0 1-0.075/100])); % saturate the upper 0.075% of the intensity range
    imshow(doseTotal', [0 0.2]); % manually define dose range of interest
    colormap(hot);
    colorbar;
    saveas(gcf, [outputDir 'DoseTotal.png']);
    
    % writeImageFormats(doseKeV, p.Reconstruction.spacing, [0 0.1], [outputDir 'DoseKeV']);
    writeMetaImage(doseKeV, [outputDir 'DoseKeV'], p.Reconstruction.spacing);
    % writeImageFormats(doseMeV, p.Reconstruction.spacing, [0 0.1], [outputDir 'DoseMeV']);
    writeMetaImage(doseMeV, [outputDir 'DoseMeV'], p.Reconstruction.spacing);
    % writeImageFormats(doseTotal, p.Reconstruction.spacing, [0 0.1], [outputDir 'DoseTotal']);
    writeMetaImage(doseTotal, [outputDir 'DoseTotal'], p.Reconstruction.spacing);
    
end

fprintf('\n');
