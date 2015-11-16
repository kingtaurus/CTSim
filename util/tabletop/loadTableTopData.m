function [sinoAtt, sinoPC, flatFieldPhotonCounts ] = loadTableTopData( dataPath, geom, flatFieldIntensity, ffieldIndexX, ffieldIndexY )
%% Load projection data from the table top system. The varian .seq
% data is firstly converted to .raw data using Viv_ind.m
% inputs:
%   dataPath - directory of the .raw projection datas
%   geom    - geometry parameters
%   flatFieldIntensity - x-ray intensity of the flat field ( default zero)
%           if zero, will measured from the reference region
%  output:
%   sinoPC - sinogram of the photon counts (intensity)
%   sinoAtt - sinogram of attenuation
%
% Meng Wu at Stanford University
% 2014 - 2015

if nargin < 3
    flatFieldIntensity = 0;
end

if nargin < 5
    ffieldIndexX = [701, 800 ];
    ffieldIndexY = [101, 500 ];
end

%% load parameteter
%dataPath = 'E:\Data\NasaFlame\June_18_Study\KrCalibration_BeamHardening_60kV20ma\Ballloon_CT_100PercentKr_60kV_20ma\';

% set the variance detector size
detWidth = 1024;
detHeight = 768;

% valid detector size defined by geometry paramters
validPixelsX = geom.detSize(1);
validPixelsY = geom.detSize(2);
noViews = geom.noViews;

fprintf('Load table top data from %s \n', dataPath);

%% show an example frame
fileName = sprintf('proj_%05i.raw', 10 );

fid=fopen([dataPath fileName] );
frame = fread(fid,[detWidth,detHeight],'float')';
fclose(fid);

validIndexX = [ detWidth/2 - validPixelsX/2 + 1 ,   detWidth/2 + validPixelsX/2];
validIndexY = [ detHeight/2 - validPixelsY/2 + 1 ,   detHeight/2 + validPixelsY/2];

if 0
    figure; imdisp( frame); hold on;
    plot( validIndexX, validIndexY );
    plot( validIndexX(2:-1:1), validIndexY );
    plot( validIndexX([1 1]), validIndexY );
    plot( validIndexX([2 2]), validIndexY );
    plot( validIndexX, validIndexY([1 1]) );
    plot( validIndexX, validIndexY([2 2]) );
    
    if flatFieldIntensity == 0
        % an air region to extrac flat field x-ray intensity
        plot( ffieldIndexX, ffieldIndexY, 'r' );
        plot( ffieldIndexX(2:-1:1), ffieldIndexY , 'r');
        plot( ffieldIndexX([1 1]), ffieldIndexY , 'r');
        plot( ffieldIndexX([2 2]), ffieldIndexY , 'r');
        plot( ffieldIndexX, ffieldIndexY([1 1]) , 'r');
        plot( ffieldIndexX, ffieldIndexY([2 2]) , 'r');
    end
end
%% load projections

sinoPC = zeros(validPixelsY, validPixelsX, noViews, 'single');
flatFieldPhotonCounts = zeros( 1, noViews );

for iv = 1:noViews
    
    if mod( iv, 50 ) == 1 || iv == 1 || iv == noViews
        fprintf('(%i/%i)... ', iv, noViews );
    end
    
    fileName = sprintf('proj_%05i.raw', iv-1 );
    
    % load the projection
    fid=fopen([dataPath fileName] );
    frame = fread(fid,[detWidth,detHeight],'float')';
    fclose(fid);
    
    sinoPC(:,:,iv) = single( frame( validIndexY(1):validIndexY(2), validIndexX(1):validIndexX(2) ) ) ;
    
    if flatFieldIntensity == 0
        % get the I0 from the reference region
        refArea = frame(  ffieldIndexY(1):ffieldIndexY(2), ffieldIndexX(1):ffieldIndexX(2) );
        flatFieldPhotonCounts(iv) = mean( refArea(:) ) ;
    end
    
end
fprintf('Done.\n');

%% compute the log of x-ray intensity to get attenuation

dummyLowPhotonCount = 10;

sinoPC( sinoPC < dummyLowPhotonCount) = dummyLowPhotonCount;
sinoAtt = zeros(validPixelsY, validPixelsX, noViews, 'single');

for iv = 1:noViews
    
    if flatFieldIntensity == 0
        flatFieldValue = flatFieldPhotonCounts(iv);
    else
        flatFieldValue = flatFieldIntensity;
    end
    
    sinoAtt(:, :, iv) = - log( sinoPC(:, :, iv) ./  flatFieldValue );
end

fprintf( 'Mean flat field intensity is %.0f. \n', mean( flatFieldValue ));

