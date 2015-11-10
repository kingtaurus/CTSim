% Read XCAT
clear all 
close all
clc

%fid = fopen('MengPhantom1_atn_1.bin'); % change filename
%phant = single( fread(fid,'single') ); % read entire phantom into a vector

%phant =  reshape(phant,281,512,601) ; % reshape


% load '1_mm_pelvis.mat';
%load 'DistortedPelvisPhantoms.mat'; % It includes Pelvis1,2,3,4 (300*700*501)
load( 'E:\MATLAB\CTData\scatter_simulation_pelvis\Pelvis5and6.mat'); % It includes Pelvis5 and 6.

% % put into standard view (patient lying on back, z = axis, z = 1 is feet)
% for i = 1:size(phant,3)
%     phant(:,:,i) = phant(:,:,i)';
% end

% truncate (if desired-- change range!)% phant = phant(40:210,80:170,:);

%in any of the images below, if you can't see the phantom, multiply phant
%by some large number (100?) before plotting

%figure; imagesc(Pelvis1(:,:,round(end/2)));

phantom = Pelvis5;
%phantom = phantom(:,:,381:480);
%phantom = phantom(:,178:510,381:480);
%phantom = upsampleVolumes( phantom, 2 );

% look at an axial slice
figure; imagesc(phantom(:,:,round(end/2)));
figure; imagesc(phantom(:,:,round(end)));
figure; imagesc(phantom(:,:,round(1)));
% look at a coronal slice
figure; imagesc(squeeze(phantom( 200,:,:)));

% look at a sagittal slice
figure;imagesc(squeeze(phantom(:,round(end/2),:)));




%%

J = single(phantom );

I = zeros( size(J), 'single' );
%I( J <= 0.0050 ) = 0;               % air
I( J <= 0.0100 ) = 1;  % lung
I( J > 0.0100 & J <= 0.0210 ) = 2;  % soft tissue
I( J > 0.0210 & J <= 0.0216 ) = 3;  % muscle
I( J > 0.0216 & J <= 0.0218 ) = 4;  % liver
I( J > 0.0218 & J <= 0.025 ) = 5;  % bone plastic
I( J > 0.025 ) = 6;                % bone cortical


K = zeros( size(J), 'single' );
%K( I == 0 ) = materialAttenuation( 60, 'air');   
K( I == 1 ) = materialAttenuation( 60, 'lung');
K( I == 2 ) = materialAttenuation( 60, 'soft_tissue'); 
K( I == 3 ) = materialAttenuation( 60, 'muscle_skeletal');   
K( I == 4 ) = materialAttenuation( 60, 'blood'); 
K( I == 5 ) = materialAttenuation( 60, 'bone_compact');
K( I == 6 ) = materialAttenuation( 60, 'bone_cortical'); 

density = J * 10 ./ K;

density( density > 1.5 ) = 1.5;

segm = I;


figure;
imagesc( squeeze( I(:,:,round(end/2)) ) );

figure;
imagesc( squeeze( density(:,:,round(end/2) )) );

%%
subRadius = 25;
subDensity = [1.1 1.05 1.03];
subTargetsDiameters = [10 8 6 4 3 2];
subTargetsRaduis = subTargetsDiameters / 2;
subDistancs = 4;

ellipsesParams = zeros(15, 6 );



for i = 1:3
    alpha = pi / 12 ;
    for j = 1:5
        ellipsesParams((i-1)*5+j,:) =  [ subRadius*cos( alpha + i*2*pi/3 ) ...
            subRadius*sin( alpha + i*2*pi/3 )+30  subTargetsRaduis(j)  subTargetsRaduis(j)    0   subDensity(i)];
        alpha = alpha + ( subTargetsRaduis(j) + subTargetsRaduis(j+1) + subDistancs ) / subRadius;
    end
    
end

ellipsesParams(:,2) = ellipsesParams(:,2) - 30;
ellipsesParams(:,1) = ellipsesParams(:,1) - 40;


meta.DimSize = size(squeeze( density(:,:,end/2))) ;
meta.ElementSpacing = [0.35 0.35];
for i =  - 10 : 10
    density(:,:,end/2+i) = addEllipses(squeeze( density(:,:,end/2+i) ), meta, ellipsesParams);
end


figure;
imagesc( squeeze( density(:,:,end/2) ) );

%%
subRadius = 12;
subDensity = [1.1 1.05 1.03];
subTargetsDiameters = [6 4 3 2];
subTargetsRaduis = subTargetsDiameters / 2;
subDistancs = 4;

ellipsesParams = zeros(12, 6 );



for i = 1:3
    alpha = pi / 12 ;
    for j = 1:3
        ellipsesParams((i-1)*3+j,:) =  [ subRadius*cos( alpha + i*2*pi/3 + pi ) ...
            subRadius*sin( alpha + i*2*pi/3 + pi )+30  subTargetsRaduis(j)  subTargetsRaduis(j)    0   subDensity(i)];
        alpha = alpha + ( subTargetsRaduis(j) + subTargetsRaduis(j+1) + subDistancs ) / subRadius;
    end
end

ellipsesParams(:,2) = ellipsesParams(:,2) - 30;
ellipsesParams(:,1) = ellipsesParams(:,1) - 40;


meta.DimSize = size(squeeze( density(:,:,end/2))) ;
meta.ElementSpacing = [0.35 0.35];
for i =  - 10 : 10
    density(:,:,end/2+i) = addEllipses(squeeze( density(:,:,end/2+i) ), meta, ellipsesParams);
end


figure;
imagesc( squeeze( density(:,:,end/2) ) );


%% save to raw data for lung imaging

S = squeeze(I(:,:,100));
D = squeeze(density(:,:,100));
meta.DimSize = size(S);
meta.ElementSpacing = [0.4 0.4];
writeMetaImage(uint8(S), 'XCATliver-median.mhd', meta);
writeMetaImage(D, 'XCATliver-median-density.mhd', meta);

meta.ElementSpacing = [0.5 0.5];
writeMetaImage(uint8(S), 'XCATliver.mhd', meta);
writeMetaImage(D, 'XCATliver-density.mhd', meta);

meta.DimSize = size(I);
meta.ElementSpacing = [0.5 0.5 0.5];
writeMetaImage(uint8(I), 'XCATliver-3d.mhd', meta);
writeMetaImage(density, 'XCATliver-3d-density.mhd', meta);


meta.DimSize = size(I);
meta.ElementSpacing = [0.4 0.4 0.4];
writeMetaImage(uint8(I), 'XCATliver-median-3d.mhd', meta);
writeMetaImage(density, 'XCATliver-median-3d-density.mhd', meta);



return;
