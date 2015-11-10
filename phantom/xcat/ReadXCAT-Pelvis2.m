% Read XCAT
clear all 
close all
clc

%fid = fopen('MengPhantom1_atn_1.bin'); % change filename
%phant = single( fread(fid,'single') ); % read entire phantom into a vector

%phant =  reshape(phant,281,512,601) ; % reshape


% load '1_mm_pelvis.mat';
% load 'DistortedPelvisPhantoms.mat'; % It includes Pelvis1,2,3,4 (300*700*501)
load( 'Pelvis7.mat'); % It includes Pelvis5 and 6.

% % put into standard view (patient lying on back, z = axis, z = 1 is feet)
% for i = 1:size(phant,3)
%     phant(:,:,i) = phant(:,:,i)';
% end
%%
% truncate (if desired-- change range!)% phant = phant(40:210,80:170,:);

%in any of the images below, if you can't see the phantom, multiply phant
%by some large number (100?) before plotting

%phantom = Pelvis1(:,:,round(end/2)); % 300*700*501 



%phantom = Pelvis1(:,:,381:480);
% phantom = Pelvis4;
phantom = Pelvis7;
%phantom = phantom(:,178:510,381:480);
%phantom = upsampleVolumes( phantom, 2 );

% look at an axial slice
figure; imagesc(phantom(:,:,round(end/2)));
figure; imagesc(phantom(:,:,round(end)));
figure; imagesc(phantom(:,:,round(1)));
% look at a coronal slice
%figure; imagesc(squeeze(phantom( 150,:,:)));
figure; imagesc(squeeze(phantom(round(end/2),:,:)))

% look at a sagittal slice
figure;imagesc(squeeze(phantom(:,round(end/2),:)));




%%

J = single(phantom);

I = zeros( size(J), 'single' );
I( J <= 0.0050 ) = 0;               % air
I( J > 0.0050 & J <= 0.0100 ) = 1;  % lung
I( J > 0.0100 & J <= 0.0250 ) = 2;  % soft tissue
%I( J > 0.0237 & J <= 0.0241 ) = 3;  % muscle        
%I( J > 0.0241 & J <= 0.0250 ) = 4;  % blood
I( J > 0.0250 & J <= 0.0300 ) = 3;  % bone plastic
I( J > 0.0300 ) = 4;                % bone cortical

K = zeros( size(J), 'single' );
K( I == 0 ) = materialAttenuation( 60, 'Air_Dry_near_sea_level');   
K( I == 1 ) = materialAttenuation( 60, 'Lung_Tissue_ICRU-44'); 
K( I == 2 ) = materialAttenuation( 60, 'Tissue_Soft_ICRU-44'); 
% K( I == 3 ) = materialAttenuation( 50, 'muscle_skeletal');   
% K( I == 4 ) = materialAttenuation( 50, 'blood'); 
K( I == 3 ) = materialAttenuation( 60, 'B-100_Bone-Equivalent_Plastic');
K( I == 4 ) = materialAttenuation( 60, 'Bone_Cortical_ICRU-44'); 

% ac_air = materialAttenuation( [1:200], 'Air_Dry_near_sea_level') / 10;   
% ac_lung = materialAttenuation( [1:200], 'Lung_Tissue_ICRU-44') / 10;   
% ac_tissue = materialAttenuation( [1:200], 'Tissue_Soft_ICRU-44') / 10;   
% ac_bone_plastic = materialAttenuation( [1:200], 'B-100_Bone-Equivalent_Plastic') / 10;   
% ac_bone_cortical = materialAttenuation( [1:200], 'Bone_Cortical_ICRU-44') / 10;  

% Relative density
density = J * 10 ./ K;

%% Display segmentation and relative density
segm = I;

figure;
imagesc( squeeze( I(:,:,round(end/2)) ) );

figure;
imagesc( squeeze( density(:,:,round(end/2 ))) );

%% save to raw data for pelvis imaging

meta.DimSize = size(I);
meta.ElementSpacing = [1.0 1.0 1.0];
writeMetaImage(uint8(I), 'XCATPelvis6-3d.mhd', meta);
writeMetaImage(density, 'XCATPelvis6-3d-density.mhd', meta);

%% Compute the true density

% read density from file
materialsDir = 'E:\MATLAB\CTSim\physicsdata\materials';

fid = fopen(fullfile(materialsDir, 'density_Air_Dry_near_sea_level.txt'), 'r');
rho = fscanf(fid, '%f');
fclose(fid);
density(I==0) = density(I==0)*rho;

fid = fopen(fullfile(materialsDir, 'density_Lung_Tissue_ICRU-44.txt'), 'r');
rho = fscanf(fid, '%f');
fclose(fid);
density(I==1) = density(I==1)*rho;

fid = fopen(fullfile(materialsDir, 'density_Tissue_Soft_ICRU-44.txt'), 'r');
rho = fscanf(fid, '%f');
fclose(fid);
density(I==2) = density(I==2)*rho;

fid = fopen(fullfile(materialsDir, 'density_B-100_Bone-Equivalent_Plastic.txt'), 'r');
rho = fscanf(fid, '%f');
fclose(fid);
density(I==3) = density(I==3)*rho;

fid = fopen(fullfile(materialsDir, 'density_Bone_Cortical_ICRU-44.txt'), 'r');
rho = fscanf(fid, '%f');
fclose(fid);
density(I==4) = density(I==4)*rho;



%% Save segmentation and density matrix for Erica
%save('segm4.mat','segm');
%save('density4.mat','density');
save('segm6.mat','segm');
save('density6.mat','density');

% %% save to raw data for lung imaging
% 
% meta.DimSize = size(I);
% meta.ElementSpacing = [1.0 1.0 1.0];
% writeMetaImage(uint8(I), 'XCATPelvis6-3d.mhd', meta);
% writeMetaImage(density, 'XCATPelvis6-3d-density.mhd', meta);


% S = squeeze(I(:,:,100));
% D = squeeze(density(:,:,100));
% meta.DimSize = size(S);
% meta.ElementSpacing = [0.4 0.4];
% writeMetaImage(uint8(S), 'XCATliver-median.mhd', meta);
% writeMetaImage(D, 'XCATliver-median-density.mhd', meta);
% 
% meta.ElementSpacing = [0.5 0.5];
% writeMetaImage(uint8(S), 'XCATliver.mhd', meta);
% writeMetaImage(D, 'XCATliver-density.mhd', meta);
% 
% meta.DimSize = size(I);
% meta.ElementSpacing = [0.5 0.5 0.5];
% writeMetaImage(uint8(I), 'XCATliver-3d.mhd', meta);
% writeMetaImage(density, 'XCATliver-3d-density.mhd', meta);
% 
% 
% meta.DimSize = size(I);
% meta.ElementSpacing = [0.4 0.4 0.4];
% writeMetaImage(uint8(I), 'XCATliver-median-3d.mhd', meta);
% writeMetaImage(density, 'XCATliver-median-3d-density.mhd', meta);



return;
