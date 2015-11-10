% Read XCAT
clear all 
close all
clc

addpath(genpath('/E/MATLAB/CTSim/phantom/xcat'))

load 'xcat_thorax_highres_60kev.mat';


phant = phant( 351:700, 101:880,  51:350 );

% look at an axial slice
figure;
imagesc(phant(:,:,end/2));

% look at a coronal slice
figure;
imagesc(squeeze(phant(end/2,:,:)));

% look at a sagittal slice
figure;
imagesc(squeeze(phant(:,end/2,:)));


%%
J = single(phant ) * 1.5 ;


% look at an axial slice
figure;
imagesc(J(:,:,end/2));


I = zeros( size(J), 'single' );
I( J <= 0.0050 ) = 0;               % air
I( J > 0.0050 & J <= 0.0100 ) = 1;  % lung
I( J > 0.0100 & J <= 0.0250 ) = 2;  % soft tissue
I( J > 0.0250 & J <= 0.0330 ) = 3;  % bone plastic
I( J > 0.0330 ) = 4;                % bone cortical


K = zeros( size(J), 'single' );
K( I == 0 ) = materialAttenuation( 60, 'Air_Dry_near_sea_level');   
K( I == 1 ) = materialAttenuation( 60, 'Lung_Tissue_ICRU-44'); 
K( I == 2 ) = materialAttenuation( 60, 'Tissue_Soft_ICRU-44'); 
K( I == 3 ) = materialAttenuation( 80, 'B-100_Bone-Equivalent_Plastic');
K( I == 4 ) = materialAttenuation( 70, 'Bone_Cortical_ICRU-44'); 

ac_air = materialAttenuation( [1:200], 'Air_Dry_near_sea_level') / 10;   
ac_lung = materialAttenuation( [1:200], 'Lung_Tissue_ICRU-44') / 10;   
ac_tissue = materialAttenuation( [1:200], 'Tissue_Soft_ICRU-44') / 10;   
ac_bone_plastic = materialAttenuation( [1:200], 'B-100_Bone-Equivalent_Plastic') / 10;   
ac_bone_cortical = materialAttenuation( [1:200], 'Bone_Cortical_ICRU-44') / 10;   

density = J * 10 ./ K;

segm = I;


figure;
imagesc( squeeze( I(:,:,end/2) ) );

figure;
imagesc( squeeze( density(:,:,end/2) ) );

return;

%% save to raw data for lung imaging

S = squeeze(I(:,:,201));
D = squeeze(density(:,:,201));
meta.DimSize = size(S);
meta.ElementSpacing = [0.6 0.6];
writeMetaImage(uint8(S), 'XCAT-thorax-median.mhd', meta);
writeMetaImage(D, 'XCAT-thorax-median-density.mhd', meta);

meta.ElementSpacing = [0.7 0.7];
writeMetaImage(uint8(S), 'XCAT-thorax.mhd', meta);
writeMetaImage(D, 'XCAT-thorax-density.mhd', meta);

meta.DimSize = size(I);
meta.ElementSpacing = [0.6 0.6 0.6];
writeMetaImage(uint8(I), 'XCAT-thorax-3d.mhd', meta);
writeMetaImage(density, 'XCAT-thorax-3d-density.mhd', meta);


meta.DimSize = size(I);
meta.ElementSpacing = [0.7 0.7 0.7];
writeMetaImage(uint8(I), 'XCAT-thorax-median-3d.mhd', meta);
writeMetaImage(density, 'XCAT-thorax-median-3d-density.mhd', meta);

