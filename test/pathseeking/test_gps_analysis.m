
%clear;
close all;

load 'liver.mat';

%% Display FBP reconstructions

ntx = 48;
nty = 90;

window0 = [0 150];
window1 = [0 150];
window2 = [-20 20];

img_fbp_sharp_slice     = img_fbp_sharp(nty+1:end-nty, ntx+1:end-ntx, ceil(end/2));
img_fbp_soft_slice      = img_fbp_soft(nty+1:end-nty, ntx+1:end-ntx, ceil(end/2));

figure; imdisp( att2hu(img_fbp_sharp_slice ), window0);
saveas( gcf, 'fbp_sharp.eps');

figure; imdisp( att2hu(img_fbp_soft_slice ), window0);
saveas( gcf, 'fbp_soft.eps');

%%
img_pwls_beta_small_slice   = img_pwls_beta_small(nty+1:end-nty, ntx+1:end-ntx, ceil(end/2));
img_pwls_beta_large_slice   = img_pwls_beta_large(nty+1:end-nty, ntx+1:end-ntx, ceil(end/2));

figure; imdisp( att2hu(img_pwls_beta_small_slice ), window0);
saveas( gcf, 'pwls_sharp.eps');

figure; imdisp( att2hu(img_pwls_beta_large_slice ), window0);
saveas( gcf, 'pwls_smooth.eps');

img_pwls_beta_mid1_slice   = img_pwls_beta_1(nty+1:end-nty, ntx+1:end-ntx, ceil(end/2));
img_pwls_beta_mid2_slice   = img_pwls_beta_2(nty+1:end-nty, ntx+1:end-ntx, ceil(end/2));

figure; imdisp( att2hu(img_pwls_beta_mid1_slice ), window0);
saveas( gcf, 'pwls_mid_1.eps');
figure; imdisp( att2hu(img_pwls_beta_mid2_slice ), window0);
saveas( gcf, 'pwls_mid_2.eps');

%%
img_mix = [img_fbp_soft_slice img_pwls_beta_small_slice img_pwls_beta_mid2_slice img_pwls_beta_large_slice];

figure; imdisp( att2hu(img_mix ), window1);
saveas( gcf, 'recon_mix.eps');


%%
close all;
indx40 = [2 9 18 24 32 40];
indx20 = [1 4 8 12 16 20];

%indx20 = [1 7 14 20];

ps_pwls  = img_seqs(nty+1:end-nty, ntx+1:end-ntx, indx20 );
ps_rog_0  = img_aps(nty+1:end-nty, ntx+1:end-ntx, indx40);
ps_rog_2  = img_tps_2(nty+1:end-nty, ntx+1:end-ntx, indx40);
ps_dog_0 = img_psadmma_1(nty+1:end-nty, ntx+1:end-ntx, indx40);
ps_dog_2 = img_psadmma_3(nty+1:end-nty, ntx+1:end-ntx, indx40);

figure; imdisp( att2hu( ps_pwls(:,:) ), window1);
text( 16, 32, 'PWLS', 'Color', 'white', 'FontSize', 24 )
saveas( gcf, 'ps_pwls.eps');

figure; imdisp( att2hu( ps_rog_0(:,:) ), window1);
text( 16, 32, 'PS-ROG-0', 'Color', 'white', 'FontSize', 24 )
saveas( gcf, 'ps_rog_0.eps');

figure; imdisp( att2hu( ps_rog_2(:,:) ), window1);
text( 16, 32, 'PS-ROG-2', 'Color', 'white', 'FontSize', 24 )
saveas( gcf, 'ps_rog_2.eps');

figure; imdisp( att2hu( ps_dog_0(:,:) ), window1);
text( 16, 32, 'PS-DOG-0', 'Color', 'white', 'FontSize', 24 )
saveas( gcf, 'ps_dog_0.eps');

figure; imdisp( att2hu( ps_dog_2(:,:) ), window1);
text( 16, 32, 'PS-DOG-2', 'Color', 'white', 'FontSize', 24 )
saveas( gcf, 'ps_dog_2.eps');


%%
% live rog 2 [3 8 17 22 25 32] dog 2 [3 9 16 24 32 40]
% thorax rog2 [1 3 12 21 27 32] dog 2 [5 8 18 25 33 40]
close all;

ps_rog_2  = img_tps_2(nty+1:end-nty, ntx+1:end-ntx, [3 8 17 22 25 32]);
ps_dog_2 = img_psadmma_3(nty+1:end-nty, ntx+1:end-ntx, [3 9 16 24 32 40]);

diff_rog_2 = att2hu(ps_pwls) - att2hu(ps_rog_2);
figure; imdisp( diff_rog_2(:,:), window2);
text( 16, 32, 'PS-ROG-2', 'Color', 'white', 'FontSize', 24 )
saveas( gcf, 'diff_rog_2.eps');


diff_dog_2 = att2hu(ps_pwls) - att2hu(ps_dog_2);
figure; imdisp( diff_dog_2(:,:), window2);
text( 16, 32, 'PS-DOG-2', 'Color', 'white', 'FontSize', 24 )
saveas( gcf, 'diff_dog_2.eps');


%%

updates_pwls = att2hu( ps_pwls(:,:,2:end)) - att2hu( ps_pwls(:,:,1:end-1) );
figure; imdisp( updates_pwls(:,:), window2);
text( 16, 32, 'PWLS', 'Color', 'white', 'FontSize', 24 )
saveas( gcf, 'updates_pwls.eps');

updates_rog_0 = att2hu(ps_rog_2(:,:,2:end)) - att2hu(ps_rog_2(:,:,1:end-1));
figure; imdisp( updates_rog_0(:,:), window2);
text( 16, 32, 'PS-ROG-0', 'Color', 'white', 'FontSize', 24 )
saveas( gcf, 'updates_rog_0.eps');

updates_rog_2 = att2hu(ps_rog_2(:,:,2:end)) - att2hu(ps_rog_2(:,:,1:end-1));
figure; imdisp( updates_rog_2(:,:), window2);
text( 16, 32, 'PS-ROG-2', 'Color', 'white', 'FontSize', 24 )
saveas( gcf, 'updates_rog_2.eps');

updates_dog_0 = att2hu(ps_dog_0(:,:,2:end)) - att2hu(ps_dog_0(:,:,1:end-1));
figure; imdisp( updates_dog_0(:,:), window2);
text( 16, 32, 'PS-DOG-0', 'Color', 'white', 'FontSize', 24 )
saveas( gcf, 'updates_dog_0.eps');

updates_dog_2 = att2hu(ps_dog_2(:,:,2:end)) - att2hu(ps_dog_2(:,:,1:end-1));
figure; imdisp( updates_dog_2(:,:), window2);
text( 16, 32, 'PS-DOG-2', 'Color', 'white', 'FontSize', 24 )
saveas( gcf, 'updates_dog_2.eps');

%%
e1 = att2hu( ps_pwls ) - att2hu(ps_dog_2);
figure; imdisp( e1 , [-20 20]);

%%
close all;

v1 = 190;
v2 = 320;
u1 = 120;
u2 = 250;

y = img_seqs(v1:v2,  u1:u2, [ 2   4    7    10   12   15   17   20] );
figure; imdisp(att2hu(y(:,:)),window1);
text( 16, 24, 'PWLS', 'FontSize', 14 )
saveas( gcf, 'roi_pwls.eps');

y = img_tps_2(v1:v2,  u1:u2, [ 1     7    12    18    23    29    34    40] );
figure; imdisp(att2hu(y(:,:)),window1);
text( 16, 24, 'PS-ROG-2', 'FontSize', 14 )
saveas( gcf, 'roi_rog_2.eps');

y = img_psadmma_3(v1:v2,  u1:u2, [ 1     7    12    18    23    29    34    40] );
figure; imdisp(att2hu(y(:,:)),window1);
text( 16, 24, 'PS-DOG-2', 'FontSize', 14 )
saveas( gcf, 'roi_dog_2.eps');

%%
v1 = nty + 1;
v2 = 360 - nty;
u1 = ntx + 1;
u2 = 512 - ntx;


d = img_pwls_beta_large_slice - img_pwls_beta_small_slice;

figure; imdisp(img_seqs(v1:v2,  u1:u2, 10),[0 1]);
type = 2;
error = zeros( 40, 20, 6 );

measurePathImages(img_pwls_beta_small_slice, img_pwls_beta_large_slice, type, d )

for j = 1 : 20
    
    t = img_seqs(v1:v2,  u1:u2, j);
    for i = 1 : 40
        
        p = img_aps(v1:v2,  u1:u2, i);
        error(i,j,1) = measurePathImages( p, t, type, d );
        
        p = img_tps_1(v1:v2,  u1:u2, i);
        error(i,j,2) = measurePathImages( p, t, type, d );
        
        p = img_tps_2(v1:v2,  u1:u2, i);
        error(i,j,3) = measurePathImages( p, t, type, d );
        
        p = img_psadmma_1(v1:v2,  u1:u2, i);
        error(i,j,4) = measurePathImages( p, t, type, d );

        p = img_psadmma_2(v1:v2,  u1:u2, i);
        error(i,j,5) = measurePathImages( p, t, type, d );
        
        p = img_psadmma_3(v1:v2,  u1:u2, i);
        error(i,j,6) = measurePathImages( p, t, type, d );
        
    end
end

%%
figure;
%plot( squeeze( error_liver_mad(:, 10 ,:) ), 'LineWidth', 2 ); hold on;
plot( squeeze( error(:, 10 ,:)), '--', 'LineWidth', 2 );
legend( 'PS-ROG-0', 'PS-ROG-1', 'PS-ROG-2', 'PS-DOG-0', 'PS-DOG-1', 'PS-DOG-2' );
xlabel('Path Seeking Frame', 'FontSize', 18 );
if type == 1
    ylabel('RMSD (HU)', 'FontSize', 18 );
else
    ylabel('MAD (HU)', 'FontSize', 18 );
end

set( gca , 'FontSize', 14);
set(gca,'units','normalized')

%%

err_min = squeeze( min( error,[], 1) );

x = 1:20;
xi = (1:40)/2;
err_min = interp1(x, err_min, xi);

figure;
%plot(  squeeze( min( error_liver_mad,[], 1) ), 'LineWidth', 2 ); hold on;
plot( squeeze( min( error,[], 1) ), '--', 'LineWidth', 2 );
legend( 'PS-ROG-0', 'PS-ROG-1', 'PS-ROG-2', 'PS-DOG-0', 'PS-DOG-1', 'PS-DOG-2' );
xlabel('PWLS Frame', 'FontSize', 18 );
if type == 1
    ylabel('RMSD (HU)', 'FontSize', 18 );
else
    ylabel('MAD (HU)', 'FontSize', 18 );
end
set( gca , 'FontSize', 14);
set(gca,'units','normalized')


%%

figure;
plot( att2hu(squeeze(x_true(246,110:130,:))'), 'LineWidth', 2 );
xlabel('Frame', 'FontSize', 20 );
ylabel('Pixel Values (HU)', 'FontSize', 20 );
set( gca , 'FontSize', 20);
axis([0 20 -50 100]);


figure;
plot( 1:2:40, att2hu(squeeze(x1_f_1_10(246,110:130,1:2:end))'), 'LineWidth', 2 );
xlabel('Frame', 'FontSize', 20 );
ylabel('Pixel Values (HU)', 'FontSize', 20 );
set( gca , 'FontSize', 20);
axis([0 40 -50 100]);