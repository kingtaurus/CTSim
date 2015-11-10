clear; clc; close all;
fprintf('Welcome to the X-ray CT simluation toolbox!\n')
fprintf('\tcopyright: Meng Wu at Stanford University (2012-2014)\n');
fprintf('\tcontact: mengwu@stanford.edu\n');

%% add path
fprintf('\tAdding needed paths... ');
addpath(genpath('../CTSim'));
addpath(genpath('../CTdata'));
fprintf('\tdone.\n');

fprintf('\tAdding output paths... ');
outputDir = [ 'output'] ;
if outputDir(end) ~= '/', outputDir = [outputDir '/']; end % make sure directory ends with a slash
if ~exist(outputDir, 'dir'), mkdir(outputDir); end % make sure output directory exists
addpath(genpath(outputDir));
fprintf('\tdone.\n');

%% Reduce dimension and downsampling

fprintf('\tRead Parameter File... ');

downsampleRate  = 2;
dimemsions      = 2;

parameters = 'parametersPCXD.ini';

p = readParametersCT(parameters);
p.Phantom.materialsFileName = [p.Phantom.materialsFileName '-3d'];

%% downsampling
if downsampleRate ~= 1
    p.Reconstruction.size       = round( p.Reconstruction.size / downsampleRate );
    p.Reconstruction.spacing    = p.Reconstruction.spacing * downsampleRate;
    p.Geometries.sizeDet    = round( p.Geometries.sizeDet / downsampleRate );
    p.Geometries.noViews    = round( p.Geometries.noViews / downsampleRate );
    p.Geometries.spacingDet = p.Geometries.spacingDet * downsampleRate;
    
end


%% reduce dimension
if dimemsions == 2
    p.Reconstruction.size       = p.Reconstruction.size(1:2);
    p.Reconstruction.spacing    = p.Reconstruction.spacing(1:2);
    p.Reconstruction.offset     = p.Reconstruction.offset(1:2);
    p.Geometries.sizeDet    = p.Geometries.sizeDet(1);
    p.Geometries.spacingDet = p.Geometries.spacingDet(1);
    p.Geometries.offsetDet  = p.Geometries.offsetDet(1);
    
    sinosDirKeV = sprintf('../CTData/SinogramDataPCXD/sinoData-%s-%s-%imm-%imm-%ix%i-%ium/', ...
        parameters(11:end-4), p.Phantom.materialsFileName, p.Geometries.SAD, p.Geometries.ADD,...
        p.Geometries.noViews, p.Geometries.sizeDet(1), round(p.Geometries.spacingDet(1)*1000));
    
else
    
    sinosDirKeV = sprintf('../CTData/SinogramDataPCXD/sinoData-%s-%s-%imm-%imm-%ix%ix%i-%ium/', ...
        parameters(11:end-4), p.Phantom.materialsFileName, p.Geometries.SAD, p.Geometries.ADD,...
        p.Geometries.noViews, p.Geometries.sizeDet(1),p.Geometries.sizeDet(2), ...
        round(p.Geometries.spacingDet(1)*1000));
    
end

save('temp.mat', 'p', 'outputDir', 'sinosDirKeV');

fprintf('\tdone.\n\n\n');


