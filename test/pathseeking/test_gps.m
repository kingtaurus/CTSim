load 'temp.mat';

% Load simulation parameters and datas
[phan, map] = loadMaterialsDensityPhantom(p);

[ geom ] = loadProjectionGeometryCT( p );

spectrum = loadSpectraCT(p, geom, 2e5);

%Compute ground trut
[imgGtAtt, imgGtHu ] = computeGroundTruth(phan, spectrum);
%figure; imdisp( imgGtAtt, [0.2 0.24] );


%% compute sinogram

[sinoRaw, tubeCurrentProfile] = simulateCTRawData( phan, geom, spectrum, sinosDir, 1, 1, 2 );

sinoAtt = processCTRawData( sinoRaw, spectrum, tubeCurrentProfile);

weights = computeWeightsPwls( sinoRaw, 0, spectrum.electronicNoise );

%% parameters
img_fbp_sharp = reconFBP( sinoAtt, geom, 'ram-lak');
imdisp( img_fbp_sharp, [0.18 0.22]   );

img_fbp_soft = reconFBP( sinoAtt, geom, 'hamming');
imdisp( img_fbp_soft, [0.18 0.22]    );

%%
nitn = 50;
numos = 20;
beta1 = 5e3;
beta2 = 2e5;
delta = 1e-3;

[img_pwls1] = reconPwlsLALMOs14( sinoAtt, weights, geom, beta1, 'huber', nitn, delta, numos );
figure; imdisp( img_pwls1,  [0.21 0.24] );

[img_pwls2] = reconPwlsLALMOs14( sinoAtt, weights, geom, beta2, 'huber', nitn, delta, numos );
figure; imdisp( img_pwls2,  [0.20 0.24]  );

[img_pwls0] = reconPwlsLALMOs14( sinoAtt, weights, geom, 2 * sqrt( beta1 * beta2 ), 'huber', nitn, delta, numos );
figure; imdisp( img_pwls0,  [0.20 0.24] );


%%
noFrames = 20;
dv = 2e-4;
p = 0.2;
os = 4;

% img_tps_0 = img_aps;
% betas_tps_0 = betas_aps;

[ img_aps, betas_aps ] = reconPwlsApproxPathSeeking( sinoAtt, weights, geom, 'huber', delta, ...
    os, numos, img_pwls_beta_small, img_pwls_beta_large, beta_small, beta_large, dv, p, noFrames, 0);
figure; imdisp( img_aps, [0.20 0.24] );

[ img_tps_0, betas_tps_0 ] = reconPwlsTruePathSeeking( sinoAtt, weights, geom, 'huber', delta, ...
    os, numos, img_pwls1, img_pwls2, beta1, beta2, dv, p, noFrames, 0);
figure; imdisp( img_tps_0, [0.20 0.24] );

[ img_tps_1, betas_tps_1 ] = reconPwlsTruePathSeeking( sinoAtt, weights, geom, 'huber', delta, ...
    os, numos, img_pwls1, img_pwls2, beta1, beta2,  dv, p, noFrames, 1);
figure; imdisp( img_tps_1, [0.20 0.24] );

[ img_tps_2, betas_tps_2 ] = reconPwlsTruePathSeeking( sinoAtt, weights, geom, 'huber', delta, ...
    os, numos, img_pwls1, img_pwls2, beta1, beta2,  dv, p, noFrames, 2);
figure; imdisp( img_tps_2, [0.20 0.24] );

%%

[ img_psadmma_1, betas_psadmm_1 ] = reconPwlsPathSeekingADMM( sinoAtt,weights, geom, 'huber', delta, ...
    numos / 2, img_pwls1, img_pwls2, beta1, beta2, noFrames, 1 );
figure; imdisp( img_psadmma_1, [0.20 0.24] );


[ img_psadmma_2, betas_psadmm_2 ] = reconPwlsPathSeekingADMM( sinoAtt,weights, geom, 'huber', delta, ...
    numos / 2, img_pwls1, img_pwls2, beta1, beta2, noFrames, 2 );
figure; imdisp( img_psadmma_2, [0.20 0.24] );


[ img_psadmma_3, betas_psadmm_3 ] = reconPwlsPathSeekingADMM( sinoAtt,weights, geom, 'huber', delta, ...
    numos / 2, img_pwls1, img_pwls2, beta1, beta2, noFrames, 3 );
figure; imdisp( img_psadmma_3, [0.20 0.24] );


%%

[ img_seqs, beta_seqs ] = reconPwlsSequencesLALM(  sinoAtt, weights, geom, [beta1 beta2], 'huber', delta, ...
    numos, nitn, noFrames / 2 );
figure; imdisp( img_seqs(:,:,:), [0.20 0.24] );

return;
%%



%imdisp(img_aps(:,:,2:2:19) - img_aps(:,:,1:2:18), [-0.002 0.002] );

%
dv = 2e-4;
p = 0.1;

img_aps_p1_dv1 = reconPwlsApproxPathSeeking( sinoAtt, weights, geom, 'huber', delta, os, img_pwls1, img_pwls2, dv, p, noFrames);
%figure; imdisp( img_aps_p1_dv1, [0.16 0.26] );

dv = 2e-4;
p = 0.2;

img_aps_p2_dv1 = reconPwlsApproxPathSeeking( sinoAtt, weights, geom, 'huber', delta, os, img_pwls1, img_pwls2, dv, p, noFrames);
%figure; imdisp( img_aps_p2_dv1, [0.16 0.26] );

dv = 4e-4;
p = 0.1;

img_aps_p1_dv2 = reconPwlsApproxPathSeeking( sinoAtt, weights, geom, 'huber', delta, os, img_pwls1, img_pwls2, dv, p, noFrames);
%figure; imdisp( img_aps_p1_dv2, [0.16 0.26] );

dv = 4e-4;
p = 0.2;

img_aps_p2_dv2 = reconPwlsApproxPathSeeking( sinoAtt, weights, geom, 'huber', delta, os, img_pwls1, img_pwls2, dv, p, noFrames);
%figure; imdisp( img_aps_p2_dv2, [0.16 0.26] );

return;
%%
figure
filename = 'aps.gif';
for n = 1 : noFrames
imdisp( att2hu( img_tps_0(:,:,n) ), [0 150] );
drawnow
frame = getframe(1);
im = frame2im(frame);
[A,map] = rgb2ind(im,256); 
	if n == 1;
		imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.2);
	else
		imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.2);
	end
end

%%
figure
filename = 'tps_diff.gif';
for n = 1 : noFrames - 2
imdisp( att2hu( img_tps_2(:,:,n+1) - img_tps_2(:,:,n)) + 1000, [-5 5] );
drawnow
frame = getframe(1);
im = frame2im(frame);
[A,map] = rgb2ind(im,256); 
	if n == 1;
		imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.3);
	else
		imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.3);
	end
end

%%


close all;

t = img_pwls_2e4;
d = img_pwls2 - img_pwls1;

k = 40;
type = 1;
error = zeros( k, 4);


measurePathImages(img_pwls2, img_pwls1, type, d )


for i = 1:40
    
    p1 = img_aps_p1_dv1(:,:,i);
    p2 = img_aps_p1_dv2(:,:,i);
    
    error(i,1) = measurePathImages(p1, t, type, d );
    error(i,2) = measurePathImages(p2, t, type, d );
    
    
    p1 = img_aps_p2_dv1(:,:,i);
    p2 = img_aps_p2_dv2(:,:,i);
    
    error(i,3) = measurePathImages(p1, t, type, d );
    error(i,4) = measurePathImages(p2, t, type, d );
    
   
end


figure;
plot(error(:,1:2:end), 'LineWidth', 2 ); hold on;
plot(error(:,2:2:end), '--', 'LineWidth', 2 ); 

legend( 'p = 10%, dv = 1 HU', 'p = 10%, dv = 2 HU','p = 20%, dv = 1 HU', 'p = 20%, dv = 2 HU' );
xlabel('Path Seeking Frame', 'FontSize', 16 );
ylabel('RMSD (HU)', 'FontSize', 16 );
set( gca , 'FontSize', 16);


%%

pwls_path_film = zeros( 256 * 10, 512 * 4 );

for j = 0 : 3
    for i = 0:9
        
        pwls_path_film( i*256+1:(i+1)*256, j*512+1:(j+1)*512 ) = img_aps( :, :, j*10+i+1);
                
    end
end

imdisp( pwls_path_film, [0.19 0.24] );