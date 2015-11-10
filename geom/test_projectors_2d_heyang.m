% Test_projectors: forwardProjectRayDriven2 (2D Matlab forward projection)
% He Yang, Nov 2014

%% Setting parameters

clear all; close all;

fprintf('\tAdding needed paths... \n');
addpath(genpath('../CTSim'));

%n= 128;
n = 500;

I = zeros(n,n);

ix = 1:n;
ix = ix + 0.5;
iix = ones(n,1) * ix;
iiy = iix';

r1 = round(n/12);
r2 = round(n/32);
r3 = round(n/3);

c1x = round(n*2/3);
c1y = round(n/2);

c2x = round(n/2);
c2y = round(n/3);

c3x = round(n/2);
c3y = round(n/2);

I( sqrt((iix-c3x).^2+(iiy-c3y).^2) < r3 ) = 1;
I( sqrt((iix-c1x).^2+(iiy-c1y).^2) < r1 ) = 2;
I( sqrt((iix-c2x).^2+(iiy-c2y).^2) < r2 ) = 5;


img = I;

figure(1)
imagesc(img); colormap gray;


geom.detSpacing = 1;
geom.detSize = 400;
%geom.detOffset = 10;
geom.detOffset = 0;


phan.size = size(img);
phan.spacing = [1.5 1.5];
%phan.offset  = [10, 10];
phan.offset = [0 0];

geom.reconSize = phan.size;
geom.reconSpacing = phan.spacing;
geom.reconOffset = phan.offset;


geom.SAD = 800;
geom.ADD = 800;
geom.SDD = geom.SAD + geom.ADD;
geom.noViews = 72;
geom.betas = (0:geom.noViews-1)' * 2*pi/geom.noViews;

%% Test projections: matlab function forwardProjectRayDriven2
tic
flatPanel = 0;
sino0 = forwardProjectRayDriven2(img, phan, geom, flatPanel);
toc

% tic
% sino2 = forwardProjectDistanceDrivenMex( img, phan, geom );
% toc
% 
% tic
% sino1 = forwardProjectRayDrivenMex( img, phan, geom );
% toc


figure(10)
imagesc(sino0); colormap gray;

flatPanel = 1;
sino1 = forwardProjectRayDriven2(img, phan, geom, flatPanel);
figure(11)
imagesc(sino1); colormap gray;
%
flatPanel = 2;
sino2 = forwardProjectRayDriven2(img, phan, geom, flatPanel);
figure(12)
imagesc(sino2); colormap gray;
% 
figure(13)
imagesc(sino0- sino2); colormap gray;
%
figure(14)
imagesc(sino1 - sino2);colormap gray;
%
figure(15)
plot(sino0(:,1)-sino2(:,1));

%% Test projectors: fbct_geom_mex3 (2d forward projection)
projType = 'proj,rd';
geom.flatPanel = 0;
tic
sino30 = fbct_geom_mex3( projType, int32(geom.reconSize), single(geom.reconSpacing), ...
     single(geom.reconOffset), int32(geom.detSize), single(geom.detSpacing),...
     single(geom.detOffset), single(geom.SAD), single(geom.ADD), single(geom.SDD),...
     int32(geom.noViews), double(geom.betas), int32(geom.flatPanel), single(img) );
toc
figure(30)
imagesc(sino30); colormap gray;
 
geom.flatPanel = 1; 
sino31 = fbct_geom_mex3( projType, int32(geom.reconSize), single(geom.reconSpacing), ...
     single(geom.reconOffset), int32(geom.detSize), single(geom.detSpacing),...
     single(geom.detOffset), single(geom.SAD), single(geom.ADD), single(geom.SDD),...
     int32(geom.noViews), double(geom.betas), int32(geom.flatPanel), single(img) );
 figure(31)
imagesc(sino31); colormap gray;
 

geom.flatPanel = 2;
sino32 = fbct_geom_mex3( projType, int32(geom.reconSize), single(geom.reconSpacing), ...
     single(geom.reconOffset), int32(geom.detSize), single(geom.detSpacing),...
     single(geom.detOffset), single(geom.SAD), single(geom.ADD), single(geom.SDD),...
     int32(geom.noViews), double(geom.betas), int32(geom.flatPanel), single(img) );
figure(32) 
imagesc(sino32); colormap gray; 

figure(33)
imagesc(sino30-sino32); colormap gray;

figure(34)
imagesc(sino31-sino32); colormap gray;

figure(35)
plot(sino30(:,1)-sino32(:,1));


%%
tic
sinoTF = forwardProjectTrapezoidicFootprintMex( img, phan, geom );
toc



figure(4)
imagesc(sinoTF); colormap gray;


figure(5)
imagesc(sino - sinoTF ); colormap gray;


tic
img = backProjectDistanceDrivenMex( sinoMex, geom );
toc

tic
imgMex = backProjectPixelDrivenMex( sinoMex, geom );
toc

figure(3)
imagesc(imgMex); colormap gray;

figure(4)
imagesc(img); colormap gray;



return;
%%








