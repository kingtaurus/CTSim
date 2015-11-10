% img = readMetaImage('./result/ballon_45/bhc/img_45kvp_air.mhd'); 

%%%read calibration phantom sinogram
  sinosDir = sprintf('../CTData/June_18_Study/KrCalibration_BeamHardening_60kV20ma/air_60kvp/sino.raw');


[ geom ] = loadProjectionGeometryCT( p );
spectrum = loadSpectraCT(p, geom, 2e6);
[sinoAtt]= loadTabletopCTSinogram( sinosDir, geom );

%%% reconstruct one slice of calibration phantom
geom.detSize = [1024 1];
% sino  = sinoAtt(135+100,:,:);
sino  = sinoAtt((end+1)/2+125,:,:);
img = reconFBP(sino,geom,'hamming');
img = convertMonoAttToHu(img,spectrum);

%%% known density and mass attenuation coefficient
density = 2.65; %gm/cc
mu_rho = XrayMu('SiO2',spectrum.energyAverage); %cm^2/gm
mu = mu_rho*density/10; %mm-1

%%% fitting func
poly2fun = 'p1*x.^2+p2*x';


theta = 0:2*pi/625:2*pi;
theta(1) = [];

geom.detSize = [1024 281];
for ii = 1:15:625 %%% downsample to save time
pathlength_map = pathlength_cal(img,theta(ii),geom);
pathlength(:,ii) = pathlength_map(round((geom.detSize(2)+1)/2-geom.detOffset(2)),:);
ii

f = fit(squeeze(sino(1,:,ii))',pathlength(:,ii)*mu,poly2fun,'Exclude',pathlength(:,ii)*mu==0)
p1(ii) = f.p1;
p2(ii) = f.p2;
end

alpha_new = [mean(p2(p2~=0)) mean(p1(p1~=0))]
save alpha_new alpha_new

%%% show BHC result
geom.detSize = [1024 1];
sino_bhc = zeros(size(sino));
for order = 1:1:length(alpha_new)
sino_bhc = sino_bhc+alpha_new(order)*sino.^(order);
end
img_bhc = reconFBP(sino_bhc,geom,'hamming');
img_bhc = convertMonoAttToHu(img_bhc,spectrum);

figure, plot(img(:,256),'r');hold on;plot(img_bhc(:,256),'b');

% save pathlength pathlength

