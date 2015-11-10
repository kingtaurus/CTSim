% Read XCAT
clear all 
close all
clc

addpath(genpath('/E/MATLAB/CTSim/phantom/xcat'))

%fid = fopen('MengPhantom1_atn_1.bin'); % change filename
%phant = single( fread(fid,'single') ); % read entire phantom into a vector

%phant =  reshape(phant,281,512,601) ; % reshape
%%


load '1_mm_updated.mat';
%load '1_mm_pelvis.mat'

% % put into standard view (patient lying on back, z = axis, z = 1 is feet)
% for i = 1:size(phant,3)
%     phant(:,:,i) = phant(:,:,i)';
% end

% truncate (if desired-- change range!)% phant = phant(40:210,80:170,:);

%in any of the images below, if you can't see the phantom, multiply phant
%by some large number (100?) before plotting

% look at an axial slice
figure;
imagesc(phantom(:,:,251));

% look at a coronal slice
figure;
imagesc(squeeze(phantom(140,:,:)));

% look at a sagittal slice
figure;
imagesc(squeeze(phantom(:,300,:)));

%%

J = single(phantom );

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

ac_air = materialAttenuation( [1:200], 'Air_Dry_near_sea_level') / 10;   
ac_lung = materialAttenuation( [1:200], 'Lung_Tissue_ICRU-44') / 10;   
ac_tissue = materialAttenuation( [1:200], 'Tissue_Soft_ICRU-44') / 10;   
ac_bone_plastic = materialAttenuation( [1:200], 'B-100_Bone-Equivalent_Plastic') / 10;   
ac_bone_cortical = materialAttenuation( [1:200], 'Bone_Cortical_ICRU-44') / 10;   

density = J * 10 ./ K;

segm = I;


figure;
imagesc( squeeze( I(:,:,200) ) );

figure;
imagesc( squeeze( density(:,:,200) ) );

return;

%% save to raw data for lung imaging

S = squeeze(I(:,:,201));
D = squeeze(density(:,:,201));
meta.DimSize = size(S);
meta.ElementSpacing = [0.8 0.8];
writeMetaImage(uint8(S), 'XCATlung-median.mhd', meta);
writeMetaImage(D, 'XCATlung-median-density.mhd', meta);

meta.ElementSpacing = [1.0 1.0];
writeMetaImage(uint8(S), 'XCATlung.mhd', meta);
writeMetaImage(D, 'XCATlung-density.mhd', meta);

meta.DimSize = size(I);
meta.ElementSpacing = [1.0 1.0 1.0];
writeMetaImage(uint8(I), 'XCATlung-3d.mhd', meta);
writeMetaImage(density, 'XCATlung-3d-density.mhd', meta);


meta.DimSize = size(I);
meta.ElementSpacing = [0.8 0.8 0.8];
writeMetaImage(uint8(I), 'XCATlung-median-3d.mhd', meta);
writeMetaImage(density, 'XCATlung-median-3d-density.mhd', meta);

%% save to raw data for shoulder imaging

S = squeeze(I(:,:,365));
D = squeeze(density(:,:,365));
meta.DimSize = size(S);

meta.ElementSpacing = [0.8 0.8];
writeMetaImage(uint8(S), 'XCATshoulder-median.mhd', meta);
writeMetaImage(D, 'XCATshoulder-median-density.mhd', meta);

meta.ElementSpacing = [1.0 1.0];
writeMetaImage(uint8(S), 'XCATshoulder.mhd', meta);
writeMetaImage(D, 'XCATshoulder-density.mhd', meta);

I1 = I(:,:,300:400);
density1 = I(:,:,300:400);

meta.DimSize = size(I1);
meta.ElementSpacing = [0.8 0.8 0.8];
writeMetaImage(uint8(I1), 'XCATshoulder-median-3d.mhd', meta);
writeMetaImage(density1, 'XCATshoulder-median-3d-density.mhd', meta);

meta.DimSize = size(I1);
meta.ElementSpacing = [1.0 1.0 1.0];
writeMetaImage(uint8(I1), 'XCATshoulder-3d.mhd', meta);
writeMetaImage(density1, 'XCATshoulder-3d-density.mhd', meta);







return;
