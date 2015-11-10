
dataPath = 'E:\Data\NeuroData\pig1\';

dataFilename = 'projdat_scan_pbv1.tif';

%projMatFilename = 'proj_matrix_moco.txt';
projMatFilename =  'projtable_01.txt';

spectraDir = '..\CTsim\physicsdata\spectra\seimens\';
spectrumName = 'spectrum_120kV.txt';

% collimation
useCollimationParameter = 0;
validPixelsX = 612;
validPixelsY = 256;
pixelWidth   = 0.616;
pixelHeight  = 0.616;
detWidth = 616;
detHeight = 480;

% recon
FOV = 215;
dummyLowPhotonCount = 10;
reconSize = [360 360 128];
reconSpacing = [0.7 0.7 1] ;

%spectrum
targetMaterial = 'W';
additionalAluminumFiltration = 21.0;






