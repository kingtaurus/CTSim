
%%% load parameter
load 'temp.mat';

dataPathAir = 'E:\Data\NasaFlame\June_18_Study\KrCalibration_BeamHardening_60kV20ma\Air_CT_60kV_20ma\';

dataPath = 'E:\Data\NasaFlame\June_18_Study\KrCalibration_BeamHardening_60kV20ma\Ballloon_CT_100PercentKr_60kV_20ma\';

geom = loadProjectionGeometryCT( p );

spectrum = loadSpectraCT(p, geom, 3156);

[~, sinoAttAir] = loadTableTopData( dataPathAir, geom );

[~, sinoAtt] = loadTableTopData( dataPath, geom );

%% Reconstruction 

% modify the reconstruction volume by 
geom.reconSize( 3 ) = 2;
geom.reconSpacing( 3 ) = 0.5;
geom.reconOffset( 3 ) = 00;

% reconstruct krypton and fuel
imgAttFBP = reconFBP( sinoAtt, geom, 'hamming' );

imgAttFBP = imgAttFBP(:,:,1);

figure, imdisp( imgAttFBP, [0 0.5] );


% reconstruct air only 
imgAttFBPAir = reconFBP( sinoAttAir, geom, 'hamming' );

imgAttFBPAir = imgAttFBPAir(:,:,1);

figure, imdisp( imgAttFBPAir, [0 0.5] );

%% Automatic rigid registration

[optimizer, metric] = imregconfig('monomodal');

[imgSubtraction, r] = imregister(imgAttFBPAir, imgAttFBP, 'rigid', optimizer, metric);

figure, imdisp( imgAttFBP - imgSubtraction, [-0.1 0.2] );


%% Two steps image registrations

% for air scan
[centersBright, radiiBright] = imfindcircles(imgAttFBPAir , [150 170], 'ObjectPolarity','bright')
T = maketform('affine', [1 0 0; 0 1 0; -(centersBright(1) - size( imgAttFBPAir, 1)/2) -(centersBright(2) - size(imgAttFBPAir, 2)/2) 1]);   %# represents translation
imgAttFBPAir = imtransform(imgAttFBPAir, T, 'XData',[1 size(imgAttFBPAir,2)], 'YData',[1 size(imgAttFBPAir,1)]);

[centersBright, radiiBright] = imfindcircles(imgAttFBPAir , [150 170], 'ObjectPolarity','bright')
% figure, imdisp( imgAir, [0 0.5] );
% viscircles(centersBright, radiiBright,'EdgeColor','b');

% for actual scan
[centersBright, radiiBright] = imfindcircles(imgAttFBP , [150 170], 'ObjectPolarity','bright')
T = maketform('affine', [1 0 0; 0 1 0; -(centersBright(1) - size( imgAttFBP, 1)/2) -(centersBright(2) - size(imgAttFBP, 2)/2) 1]);   %# represents translation
imgAttFBP = imtransform(imgAttFBP, T, 'XData',[1 size(imgAttFBP,2)], 'YData',[1 size(imgAttFBP,1)]);

% [centersBright, radiiBright] = imfindcircles(imgAttFBP , [150 170], 'ObjectPolarity','bright');
% figure, imdisp( imgAttFBPAir, [0 0.5] );
% viscircles(centersBright, radiiBright,'EdgeColor','b');

rotAngles = -10:0.2:0;
correctionScores = zeros( 1, length( rotAngles ));

% find the optimal re
for i = 1 : length( rotAngles )
    img = imrotate( imgAttFBPAir, rotAngles(i)/180*pi, 'bilinear', 'crop');
    correctionScores(i) = norm( imgAttFBP - img, 2);
end

figure;
plot( rotAngles, correctionScores );

[~, j] = min( correctionScores );
img = imrotate( imgAttFBPAir, rotAngles(j)/180*pi, 'bilinear', 'crop');
figure, imdisp( imgAttFBP - img, [-0.5 0.5] );


%% The following are Tony's original codes 2015.
% I changed mainly for more efficient immplementation. 

img_Kr = readMetaImage('./result/ballon_45/img_45kvp_100Kr_burner.mhd');
img_air = readMetaImage('./result/ballon_45/img_45kvp_air.mhd');

img_Kr_central = img_Kr(:,:,end/2);
img_air_central = img_air(:,:,end/2);

img_Kr_central = img_Kr_central.*(img_Kr_central>0);
img_air_central = img_air_central.*(img_air_central>0);

roi_x = 175:340;
roi_y = 160:335;

theta = -10:0.01:10;
norm_num = 2;

for ii = 1:1:length(theta)
    img_Kr_central_rt = imrotate(img_Kr_central,theta(ii)/180*pi,'bicubic','crop');
    norm_value(ii) = norm(img_Kr_central_rt(roi_x,roi_y)-img_air_central(roi_x,roi_y),norm_num);
end

[min_norm_value, min_idx] = min(norm_value);
theta(min_idx)
img = imrotate(img_Kr,theta(min_idx)/180*pi,'bilinear','crop');
 writeMetaImage(img,'./result/ballon_45/img_45kvp_100Kr_burner_rot.mhd');
