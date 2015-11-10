clear; clc; close all;
fprintf('Welcome to the X-ray CT simluation toolbox!\n')
fprintf('\tcopyright: Meng Wu at Stanford University (2012-2014)\n');
fprintf('\tcontact: mengwu@stanford.edu\n');

%% add path
fprintf('\tAdding needed paths... ');
addpath(genpath('../CTSim'));
addpath(genpath('../CTData'));
%addpath(genpath('../Tools'));
fprintf('\tdone.\n');

fprintf('\tAdding output paths... ');
outputDir = [ 'output'] ;
if outputDir(end) ~= '/', outputDir = [outputDir '/']; end % make sure directory ends with a slash
if ~exist(outputDir, 'dir'), mkdir(outputDir); end % make sure output directory exists
addpath(genpath(outputDir));
fprintf('\tdone.\n');

%% Reduce dimension and downsampling

downsampleRate  = 1;
dimemsions      = 3;

 parameters = 'parametersTabletop.ini';
% parameters = 'parametersTrueBeamTomo.ini';
% parameters = 'parametersBeamHardening.ini';
% parameters = 'parametersSiemensZeego.ini';
% parameters = 'parametersGElightspeed.ini';
% parameters = 'parametersVarianCBCT.ini';
% parameters = 'parametersSOMATOM.ini'; 

fprintf('\tRead Parameter File... ');

p = readParametersCT(parameters);


%downsampling
if downsampleRate ~= 1
    p.Reconstruction.size       = round( p.Reconstruction.size / downsampleRate );
    p.Reconstruction.spacing    = p.Reconstruction.spacing * downsampleRate;
    p.Geometries.sizeDet        = round( p.Geometries.sizeDet / downsampleRate );
    p.Geometries.noViews        = round( p.Geometries.noViews / downsampleRate );
    p.Geometries.spacingDet     = p.Geometries.spacingDet * downsampleRate;
    
end

%reduce dimension
if dimemsions == 2
    p.Reconstruction.size       = p.Reconstruction.size(1:2);
    p.Reconstruction.spacing    = p.Reconstruction.spacing(1:2);
    p.Reconstruction.offset     = p.Reconstruction.offset(1:2);
    p.Geometries.sizeDet    = p.Geometries.sizeDet(1);
    p.Geometries.spacingDet = p.Geometries.spacingDet(1);
    p.Geometries.offsetDet  = p.Geometries.offsetDet(1);
    
    sinosDir = sprintf('../CTData/SinogramDataCT/sinoData-%s-%s-%imm-%imm-%ix%i-%ium/', ...
        parameters(11:end-4), p.Phantom.materialsFileName, p.Geometries.SAD, p.Geometries.ADD,...
        p.Geometries.noViews, p.Geometries.sizeDet(1), round(p.Geometries.spacingDet(1)*1000));
    
else
    p.Phantom.materialsFileName = [p.Phantom.materialsFileName '-3d'];
    sinosDir = sprintf('../CTData/SinogramDataCT/sinoData-%s-%s-%imm-%imm-%ix%ix%i-%ium/', ...
        parameters(11:end-4), p.Phantom.materialsFileName, p.Geometries.SAD, p.Geometries.ADD,...
        p.Geometries.noViews, p.Geometries.sizeDet(1),p.Geometries.sizeDet(2), ...
        round(p.Geometries.spacingDet(1)*1000));
    
end

save('temp.mat', 'p', 'outputDir', 'sinosDir');

fprintf('\tdone.\n\n\n');


