clear all
close all 
clc

filename = 'Meng_alignment_image.dat_0010'; % modify image file name
fid = fopen(filename);
tline = fgets(fid);
xpix = 641;% 642
zpix = 64; % 64

% run through and discard the header
for i = 1:18 % this is 17 for a single projection, 18 if one of a series of images
       tline = fgets(fid);
end

primary = zeros(xpix,zpix);
compton = zeros(xpix,zpix);
rayleigh = zeros(xpix,zpix);
multiscatter = zeros(xpix,zpix);
curved_primary = zeros(xpix,zpix);
curved_all = zeros(xpix,zpix);
curved_grid = zeros(xpix,zpix);
bin_5keV = zeros(xpix,zpix);
bin_10keV = zeros(xpix,zpix);
bin_15keV = zeros(xpix,zpix);
bin_20keV = zeros(xpix,zpix);
bin_25keV = zeros(xpix,zpix);
bin_30keV = zeros(xpix,zpix);
bin_35keV = zeros(xpix,zpix);
bin_40keV = zeros(xpix,zpix);
bin_45keV = zeros(xpix,zpix);
bin_50keV = zeros(xpix,zpix);
bin_55keV = zeros(xpix,zpix);
bin_60keV = zeros(xpix,zpix);
bin_65keV = zeros(xpix,zpix);
bin_70keV = zeros(xpix,zpix);
bin_75keV = zeros(xpix,zpix);
bin_80keV = zeros(xpix,zpix);
bin_85keV = zeros(xpix,zpix);
bin_90keV = zeros(xpix,zpix);
bin_95keV = zeros(xpix,zpix);
bin_100keV = zeros(xpix,zpix);
bin_105keV = zeros(xpix,zpix);
bin_110keV = zeros(xpix,zpix);
bin_115keV = zeros(xpix,zpix);
bin_120keV = zeros(xpix,zpix);

for j = 1:zpix
    for i = 1:xpix
        row = fscanf(fid, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f \n', 31);
        %row = fscanf(fid, '%f %f %f %f %f %f %f \n', 7);
        %row = fscanf(fid, '%f %f %f %f %f \n', 4);
        primary(i,j) = row(1);
        compton(i,j) = row(2);
        rayleigh(i,j) = row(3);
        multiscatter(i,j) = row(4);
        curved_primary(i,j) = row(5);
        curved_all(i,j) = row(6);
        curved_grid(i,j) = row(7);
        bin_5keV(i,j) = row(8);
        bin_10keV(i,j) = row(9);
        bin_15keV(i,j) = row(10);
        bin_20keV(i,j) = row(11);
        bin_25keV(i,j) = row(12);
        bin_30keV(i,j) = row(13);
        bin_35keV(i,j) = row(14);
        bin_40keV(i,j) = row(15);
        bin_45keV(i,j) = row(16);
        bin_50keV(i,j) = row(17);
        bin_55keV(i,j) = row(18);
        bin_60keV(i,j) = row(19);
        bin_65keV(i,j) = row(20);
        bin_70keV(i,j) = row(21);
        bin_75keV(i,j) = row(22);
        bin_80keV(i,j) = row(23);
        bin_85keV(i,j) = row(24);
        bin_90keV(i,j) = row(25);
        bin_95keV(i,j) = row(26);
        bin_100keV(i,j) = row(27);
        bin_105keV(i,j) = row(28);
        bin_110keV(i,j) = row(29);
        bin_115keV(i,j) = row(30);
        bin_120keV(i,j) = row(31);
    end
    row = fscanf(fid, '\n',0);
end

% Look at original result of MC-GPU on flat panel detector

% figure
% imagesc(primary)
% caxis([0,7])
% title('primary')
% colorbar

% Look at primary image on curved detector

figure
imagesc(curved_primary)
caxis([0,7])
title('curved primary')
%colorbar

% Look at image with primary and scatter on curved detector

figure
imagesc(curved_all)
caxis([0,7])
title('curved all')

curved_grid(curved_grid > 0) = 3;
curved_grid(curved_grid < 0) = 0;

figure
imagesc(curved_grid)
caxis([0,7])
title('curved grid')
% 
% figure
% imagesc(bin_5keV+bin_10keV+bin_15keV+bin_20keV+bin_25keV+bin_30keV+bin_35keV+bin_40keV+bin_45keV+bin_50keV+bin_55keV+bin_60keV+bin_65keV...
%     +bin_70keV+bin_75keV+bin_80keV+bin_85keV+bin_90keV+bin_95keV+bin_100keV+bin_105keV+bin_110keV+bin_115keV+bin_120keV+bin_5keV)
% caxis([0,7])
% title('summed bins')
% 
% figure
% imagesc(bin_20keV+bin_15keV+bin_10keV+bin_5keV)
% title('0-20 keV')
% caxis([0 7])
% 
% figure
% imagesc(bin_40keV+bin_35keV+bin_30keV+bin_25keV)
% title('20-40 keV')
% caxis([0 7])
% 
% figure
% imagesc(bin_60keV+bin_55keV+bin_50keV+bin_45keV)
% title('40-60 keV')
% caxis([0 7])
% 
% figure
% imagesc(bin_80keV+bin_75keV+bin_70keV+bin_65keV)
% title('60-80 keV')
% caxis([0 7])
% 
% figure
% imagesc(bin_100keV+bin_95keV+bin_90keV+bin_85keV)
% title('80-100 keV')
% caxis([0 7])
% 
% figure
% imagesc(bin_120keV+bin_115keV+bin_110keV+bin_105keV)
% title('100-120 keV')
% caxis([0 7])
% 
% 
curved_SPR = (curved_all-curved_primary)./curved_primary;
curved_grid_SPR = (curved_grid-curved_primary)./curved_primary;

figure
imagesc(curved_SPR)
caxis([0,2])
title('curved SPR')

% threshold
%curved_grid_SPR(curved_grid_SPR >= 0) = 1;
%curved_grid_SPR(curved_grid_SPR <= 0) = 0;

figure
imagesc(curved_grid_SPR)
caxis([0,2])
title('curved grid SPR')


    