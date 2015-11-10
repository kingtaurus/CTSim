function [curved_primary, curved_all, curved_grid] = readScatterSimulationResults( filename, dims )


fid = fopen(filename);
tline = fgets(fid);
xpix = dims(2);
zpix = dims(1);

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

end