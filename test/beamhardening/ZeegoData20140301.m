
dataPath = 'E:\MATLAB\CTData\zeegoHeadScan\';

dataFilename = 'projection_2.IMA';

projMatFilename = 'pmatix\projtable_01.txt';

tubeCurretFilename = 'current_2.txt';

exposureTimeFilename = 'time_2.txt';


spectraDir = '..\CTsim\physicsdata\spectra\seimens\';
spectrumName = 'spectrum_90kV.txt';

% collimation
useCollimationParameter = 0;
validPixelsX = 1200;
validPixelsY = 320;
pixelWidth   = 0.308;
pixelHeight  = 0.308;

% recon
FOV = 250;
dummyLowPhotonCount = 10;
reconSize = [512 512 64];
reconSpacing = [0.4 0.4 1];

%spectrum
targetMaterial = 'W';
additionalAluminumFiltration = 21.0;






