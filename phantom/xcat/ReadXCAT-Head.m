% Read XCAT
clear all
close all
clc

fid = fopen('MengPhantom1_atn_1.bin'); % change filename
phant = single( fread(fid,'single') ); % read entire phantom into a vector

%%
phant =  reshape(phant,512,512,601) ; % reshape


% put into standard view (patient lying on back, z = axis, z = 1 is feet)
for i = 1:size(phant,3)
    phant(:,:,i) = phant(:,:,i)';
end

% truncate (if desired-- change range!)% phant = phant(40:210,80:170,:);

%in any of the images below, if you can't see the phantom, multiply phant
%by some large number (100?) before plotting

% look at an axial slice
imagesc(phant(:,:,500))

% look at a coronal slice
slice = reshape(phant(300,:,:),512,601,1);
imagesc(slice(:,:,1));

% look at a sagittal slice
%slice = reshape(phant(:,120,:),512,601,1);
%imagesc(slice(:,:,1));

%%

head = zeros([ 201,151,141], 'single' );

head(:,:,1:end) = single( phant(150:350,180:330,461:601));

J = upsampleVolumes( head, 3 );


figure;
imagesc( squeeze( J(:,:,end/2)-1 ) );

%%

I = zeros( size(J), 'single' );
I( J <= 0.0237 ) = 1;  % soft tissue
I( J > 0.0237 & J <= 0.02391 ) = 2;  % brain
I( J > 0.0239 & J <= 0.02501 ) = 3;  % muscle
I( J > 0.02501 ) = 4;                % bone cortical


K = zeros( size(J), 'single' );
K( I == 1 ) = materialAttenuation( 50, 'soft_tissue');
K( I == 2 ) = materialAttenuation( 50, 'brain');
K( I == 3 ) = materialAttenuation( 50, 'muscle_striated');
K( I == 4 ) = materialAttenuation( 80, 'bone_cortical');


density = J * 10 ./ K;


figure;
imagesc( squeeze( I(:,:,end/2) ) );

density( density > 0.9 & density < 1.1 ) = 1;
figure;
imagesc( squeeze( density(:,:,end/2) ) );


figure;
imagesc( squeeze( I(end/2,:,:) ) );

return;

%%
subRadius = 18;
subDensity = [1.1 1.05 1.03];
subTargetsDiameters = [8 6 4 3 2];
subTargetsRaduis = subTargetsDiameters / 2;
subDistancs = 4;

ellipsesParams = zeros(12, 6 );



for i = 1:3
    alpha = pi / 12 ;
    for j = 1:4
        ellipsesParams((i-1)*4+j,:) =  [ subRadius*cos( alpha + i*2*pi/3 ) ...
            subRadius*sin( alpha + i*2*pi/3 )+30  subTargetsRaduis(j)  subTargetsRaduis(j)    0   subDensity(i)];
        alpha = alpha + ( subTargetsRaduis(j) + subTargetsRaduis(j+1) + subDistancs ) / subRadius;
    end
    
end

ellipsesParams(:,2) = ellipsesParams(:,2) - 5;

meta.DimSize = size(squeeze( density(:,:,end/2))) ;
meta.ElementSpacing = [0.35 0.35];
for i =  - 10 : 10
    density(:,:,end/2+i) = addEllipses(squeeze( density(:,:,end/2+i) ), meta, ellipsesParams);
end


figure;
imagesc( squeeze( density(:,:,end/2) ) );

%% save to raw data
meta.DimSize = size(J);
meta.ElementSpacing = [0.35 0.35 0.35];
writeMetaImage(uint8(I), 'XCAThead-3d.mhd', meta);
writeMetaImage(density, 'XCAThead-3d-density.mhd', meta);
