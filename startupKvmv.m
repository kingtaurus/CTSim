%% Initial Environment for KV/MV Project

%% Clear MATLAB

clear; clc; close all;
fprintf('Welcome to the X-ray CT simluation toolbox!\n')
fprintf('\tcopyright: Meng Wu at Stanford University (2012-2014)\n');
fprintf('\tcontact: mengwu@stanford.edu\n');

%% Add path
fprintf('\tAdding needed paths... ');
addpath(genpath('../CTSim'));
addpath(genpath('../CTData'));
fprintf('\tdone.\n');

fprintf('\tAdding output paths... ');
outputDir = [ 'output'] ;
if outputDir(end) ~= '/', outputDir = [outputDir '/']; end % make sure directory ends with a slash
if ~exist(outputDir, 'dir'), mkdir(outputDir); end % make sure output directory exists
addpath(genpath(outputDir));
fprintf('\tdone.\n');

%% Load system parameters

downsampleRate = 1;
dimemsions = 3;

% parameters = 'parametersDentalKvmv.ini';
% parameters = 'parametersHeadKvmv.ini';
% parameters = 'parametersHipKvmv.ini';
 parameters = 'parametersTrueBeam.ini';

fprintf('\tRead Parameter File... ');
p = readParametersKVMV(parameters);
p.Phantom.materialsFileName = [p.Phantom.materialsFileName '-3d'];

%% Downsampling

if downsampleRate ~= 1
    p.Reconstruction.size       = round( p.Reconstruction.size / downsampleRate );
    p.Reconstruction.spacing    = p.Reconstruction.spacing * downsampleRate;
    p.Geometries.keV.sizeDet    = round( p.Geometries.keV.sizeDet / downsampleRate );
    p.Geometries.keV.spacingDet = p.Geometries.keV.spacingDet * downsampleRate;
    p.Geometries.MeV.sizeDet    = round( p.Geometries.MeV.sizeDet / downsampleRate );
    p.Geometries.MeV.spacingDet = p.Geometries.MeV.spacingDet * downsampleRate;
end

%% Reduce dimension

if dimemsions == 2
    p.Reconstruction.size       = p.Reconstruction.size(1:2);
    p.Reconstruction.spacing    = p.Reconstruction.spacing(1:2);
    p.Reconstruction.offset     = p.Reconstruction.offset(1:2);
    p.Geometries.keV.sizeDet    = p.Geometries.keV.sizeDet(1);
    p.Geometries.keV.spacingDet = p.Geometries.keV.spacingDet(1);
    p.Geometries.keV.offsetDet  = p.Geometries.keV.offsetDet(1);
    p.Geometries.MeV.sizeDet    = p.Geometries.MeV.sizeDet(1);
    p.Geometries.MeV.spacingDet = p.Geometries.MeV.spacingDet(1);
    p.Geometries.MeV.offsetDet  = p.Geometries.MeV.offsetDet(1);
    
    sinosDirKeV = sprintf('../CTData/SinogramDataKvmv/sinoData-%s-%imm-%imm-%ix%i-%ium/', ...
        p.Phantom.materialsFileName, p.Geometries.keV.SAD, p.Geometries.keV.ADD,...
        p.Geometries.keV.noViews, p.Geometries.keV.sizeDet(1),...
        round(p.Geometries.keV.spacingDet(1)*1000));
    
    sinosDirMeV = sprintf('../CTData/SinogramDataKvmv/sinoData-%s-%imm-%imm-%ix%i-%ium/', ...
        p.Phantom.materialsFileName, p.Geometries.MeV.SAD, p.Geometries.MeV.ADD,...
        p.Geometries.MeV.noViews, p.Geometries.MeV.sizeDet(1),...
        round(p.Geometries.MeV.spacingDet(1)*1000));
    
else
    
    sinosDirKeV = sprintf('../CTData/SinogramDataKvmv/sinoData-%s-%imm-%imm-%ix%ix%i-%ium/', ...
        p.Phantom.materialsFileName, p.Geometries.keV.SAD, p.Geometries.keV.ADD,...
        p.Geometries.keV.noViews, p.Geometries.keV.sizeDet(1),p.Geometries.keV.sizeDet(2),...
        round(p.Geometries.keV.spacingDet(1)*1000));
    
    sinosDirMeV = sprintf('../CTData/SinogramDataKvmv/sinoData-%s-%imm-%imm-%ix%ix%i-%ium/', ...
        p.Phantom.materialsFileName, p.Geometries.MeV.SAD, p.Geometries.MeV.ADD,...
        p.Geometries.MeV.noViews, p.Geometries.MeV.sizeDet(1),p.Geometries.MeV.sizeDet(2),...
        round(p.Geometries.MeV.spacingDet(1)*1000));
end

save('temp.mat', 'p', 'outputDir', 'sinosDirKeV', 'sinosDirMeV');

fprintf('\tdone.\n\t');

%% Initialize varian's file reader
addpath(genpath('../Tools/IALMLToolbox'));
addpath(genpath('../Tools/XIMToolbox'));
startupIALML;



