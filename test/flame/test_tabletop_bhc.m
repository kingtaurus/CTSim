
%%% load parameter
load 'temp.mat';


[ geom ] = loadProjectionGeometryCT( p );
spectrum = loadSpectraCT(p, geom, 2e6);
[sinoAtt]= loadTabletopCTSinogram( sinosDir, geom );

% geom.betas = geom.betas-0.1335;
% sinoAtt = circshift(sinoAtt,[0 0 -13]);
%% Normal X-ray CT


%% reconstruction

imgAttFBP = reconFBP( sinoAtt, geom, 'hamming' );
imgFBP = convertMonoAttToHu( imgAttFBP, spectrum) ;


figure, imshow( imgAttFBP(:,:), [0 0.6] ); %colormap gray;
% writeMetaImage(imgAttFBP,'./result/flame/KryptonJetNoTubePaper_avg.mhd');
writeMetaImage(imgFBP,'./result/ballon_60/test.mhd');

%%% Beam hardening correction
%%% An iterative approach to the beam hardening correction in cone beam CT
%%% Med Phys. 2000 Jan;27(1):23-9. Hsieh J, et al.

if 1
%%% pre-correction for quartz tube
temp = load('alpha_new_quartz.mat');
alpha = temp.alpha_new;
sinoAtt_bhc = zeros(size(sinoAtt));
for order = 1:1:length(alpha)
sinoAtt_bhc = sinoAtt_bhc+alpha(order)*sinoAtt.^(order);
end

%%% correction for porous media
clear sino_porous map_porous;
temp = load('alpha_new_porous.mat');
beta = temp.alpha_new;
%%%segmentation
map_porous(:,:,:) = imgAttFBP.*(imgFBP>0&imgFBP<1000);
theta = 0:2*pi/625:2*pi;
theta(1) = [];
detOffset = geom.detOffset;
% geom.detOffset = [geom.detOffset(1) 0];
%%% forward-projection
sino_porous = forwardProjectMex(map_porous,geom);
% for ii = 1:1:625
% sino_porous(:,:,ii) = forward_proj(map_porous,theta(ii),geom)/10;
% end
% sino_porous = reshape(sino_porous,[1 1024 625]);
delta_sino = (alpha(1)-beta(1))*sino_porous+(alpha(2)-beta(2))*sino_porous.^2;
sinoAtt_bhc = sinoAtt_bhc-delta_sino;

geom.detOffset = detOffset;
imgAttFBP_bhc = reconFBP(sinoAtt_bhc,geom,'hamming');
imgFBP_bhc = convertMonoAttToHu(imgAttFBP_bhc,spectrum);
writeMetaImage(imgFBP_bhc,'./result/ballon_60/bhc/test_bhc.mhd')

figure, plot(imgFBP(:,256),'r');hold on;plot(imgFBP_bhc(:,256),'b');

end

return;