
dataPath = 'E:\Data\NeuroData\patient1\';

dataFilename = 'proj_dat_patient1_fill1_1.tif';

%projMatFilename = 'proj_matrix_moco.txt';
projMatFilename =  'projtable_original.txt';

spectraDir = '..\CTsim\physicsdata\spectra\seimens\';
spectrumName = 'spectrum_90kV.txt';

% collimation
useCollimationParameter = 0;
validPixelsX = 612;
validPixelsY = 256;
pixelWidth   = 0.616;
pixelHeight  = 0.616;
detWidth = 616;
detHeight = 480;

% recon
FOV = 220;
dummyLowPhotonCount = 10;
reconSize = [360 360 128];
reconSpacing = [0.65 0.65 0.9] ;

%spectrum
targetMaterial = 'W';
additionalAluminumFiltration = 21.0;






