% Compute the average first

numTest = 7;
numRot  = 4;

sliceIdx = [9 12 15 18 21];
numSlice = length( sliceIdx );

%% first get the average images

disp('Compute the average images for FBP and PWLS reconstruction.');

img_fbp_soft_avg = 0;
img_fbp_sharp_avg = 0;

img_pwls_avg_1 = 0;
img_pwls_avg_2 = 0;
img_pwls_avg_3 = 0;
img_pwls_avg_4 = 0;

for i = 1 : numTest
    % load the image
    load(['test_100mas_' num2str(i) '.mat' ]);
    % roate the image by 0, 90, 180, and 270 degrees to increase sampling rate
    for k = 0 : numRot - 1
        % fbp reconstruction
        img_fbp_soft_avg    = img_fbp_soft_avg + rot90( img_fbp_soft(:,:,end/2), k);
        img_fbp_sharp_avg   = img_fbp_sharp_avg + rot90( img_fbp_sharp(:,:,end/2), k);
        % pwls reconstruction
        img_pwls_avg_1   = img_pwls_avg_1 + rot90( img_pwls_beta_small(:,:,end/2), k);
        img_pwls_avg_2   = img_pwls_avg_2 + rot90( img_pwls_beta_1(:,:,end/2), k);
        img_pwls_avg_3   = img_pwls_avg_3 + rot90( img_pwls_beta_2(:,:,end/2), k);
        img_pwls_avg_4   = img_pwls_avg_4 + rot90( img_pwls_beta_large(:,:,end/2), k);
    end
end

img_fbp_soft_avg   = img_fbp_soft_avg  / (numTest * numRot );
img_fbp_sharp_avg  = img_fbp_sharp_avg / (numTest * numRot );

img_pwls_avg_1 =  img_pwls_avg_1  / (numTest * numRot );
img_pwls_avg_2 =  img_pwls_avg_2  / (numTest * numRot );
img_pwls_avg_3 =  img_pwls_avg_3  / (numTest * numRot );
img_pwls_avg_4 =  img_pwls_avg_4  / (numTest * numRot );

%% now compute the noise power spectrum

disp('Compute the noise power spectra for FBP and PWLS reconstruction.');

img_fbp_soft_nps = 0;
img_fbp_sharp_nps = 0;

img_pwls_nps_1 = 0;
img_pwls_nps_2 = 0;
img_pwls_nps_3 = 0;
img_pwls_nps_4 = 0;

for i = 1 : numTest
    % load the image
    load(['test_100mas_' num2str(i) '.mat' ]);
    
    % roate the image by 0, 90, 180, and 270 degrees to increase sampling rate
    for k = 0 : numRot - 1
        % fbp reconstruction
        img_diff = rot90( img_fbp_soft(:,:,end/2), k) - img_fbp_soft_avg;
        img_fbp_soft_nps = img_fbp_soft_nps + fftshift( abs( fft2( img_diff( 129: 384, 129:384)) ).^2 ) ;
        
        img_diff = rot90( img_fbp_sharp(:,:,end/2), k) - img_fbp_sharp_avg;
        img_fbp_sharp_nps = img_fbp_sharp_nps + fftshift( abs( fft2( img_diff( 129: 384, 129:384)) ).^2 ) ;
        
        % pwls reconstruction
        img_diff = rot90( img_pwls_beta_small(:,:,end/2), k) - img_pwls_avg_1;
        img_pwls_nps_1 = img_pwls_nps_1 + fftshift( abs( fft2( img_diff( 129: 384, 129:384)) ).^2 ) ;
        
        img_diff = rot90( img_pwls_beta_1(:,:,end/2), k) - img_pwls_avg_2;
        img_pwls_nps_2 = img_pwls_nps_2 + fftshift( abs( fft2( img_diff( 129: 384, 129:384)) ).^2 ) ;
        
        img_diff = rot90( img_pwls_beta_2(:,:,end/2), k) - img_pwls_avg_3;
        img_pwls_nps_3 = img_pwls_nps_3 + fftshift( abs( fft2( img_diff( 129: 384, 129:384)) ).^2 ) ;
        
        img_diff = rot90( img_pwls_beta_large(:,:,end/2), k) - img_pwls_avg_4;
        img_pwls_nps_4 = img_pwls_nps_4 + fftshift( abs( fft2( img_diff( 129: 384, 129:384)) ).^2 ) ;
        
    end
end

img_fbp_soft_nps = img_fbp_soft_nps / (numTest * numRot);
img_fbp_sharp_nps = img_fbp_sharp_nps / (numTest * numRot);
img_pwls_nps_1 = img_pwls_nps_1 / (numTest * numRot);
img_pwls_nps_2 = img_pwls_nps_2 / (numTest * numRot);
img_pwls_nps_3 = img_pwls_nps_3 / (numTest * numRot);
img_pwls_nps_4 = img_pwls_nps_4 / (numTest * numRot);

%% Show FBP image and NPS
disp('Show NPS for FBP and PWLS reconstruction.');

figure; imdisp( img_fbp_soft_nps );
figure; imdisp( img_fbp_sharp_nps );

% Show PWLS image and NPS
figure; imdisp( img_pwls_nps_1 );
figure; imdisp( img_pwls_nps_2 );
figure; imdisp( img_pwls_nps_3 );
figure; imdisp( img_pwls_nps_4 );

%% now for the path images
disp('Compute the average images for path seeking.');

img_path_avg_1 = zeros( 512, 512, 40 );
img_path_avg_2 = zeros( 512, 512, 40 );
for i = 1 : numTest
    % load the image
    load(['test_100mas_' num2str(i) '.mat' ]);
    for path_id = 1:40      
        % roate the image by 0, 90, 180, and 270 degrees to increase sampling rate
        for k = 0 : numRot - 1
            img_path_avg_1(:,:,path_id)   = img_path_avg_1(:,:,path_id) + rot90( img_ps1(:,:,path_id), k);
            img_path_avg_2(:,:,path_id)   = img_path_avg_2(:,:,path_id) + rot90( img_ps2(:,:,path_id), k);
        end
    end
end
img_path_avg_1 = img_path_avg_1 /  (numTest * numRot);
img_path_avg_2 = img_path_avg_2 /  (numTest * numRot);


%%

disp('Compute the NPS for path seeking.');

img_path_nps_1 = zeros( 256, 256, 40 );
img_path_nps_2 = zeros( 256, 256, 40 );
for i = 1 : numTest
    % load the image
    load(['test_100mas_' num2str(i) '.mat' ]);
    for path_id = 1:40      
        % roate the image by 0, 90, 180, and 270 degrees to increase sampling rate
        for k = 0 : numRot - 1
            
            img_diff = rot90( img_ps1(:,:,path_id), k) - img_path_avg_1(:,:,path_id);
            img_path_nps_1(:,:,path_id) = img_path_nps_1(:,:,path_id) + fftshift( abs( fft2( img_diff( 129:384, 129:384)) ).^2 ) ;
            
            img_diff = rot90( img_ps2(:,:,path_id), k) - img_path_avg_2(:,:,path_id);
            img_path_nps_2(:,:,path_id) = img_path_nps_2(:,:,path_id) + fftshift( abs( fft2( img_diff( 129:384, 129:384)) ).^2 ) ;
        end
    end
end
img_path_nps_1 = img_path_nps_1 /  (numTest * numRot) * ( 0.7^2 / 256^2 );
img_path_nps_2 = img_path_nps_2 /  (numTest * numRot) * ( 0.7^2 / 256^2 );

%%
close all;

img_path_nps_all_norm = img_path_nps_2;
for i = 1:40
    img_path_nps_all_norm( :,:, i ) = img_path_nps_all_norm(:,:,i) / max( max( img_path_nps_all_norm(:,:,i) ) );    
end

img_path_nps_all_norm = img_path_nps_all_norm(:,:,1:4:30);
figure; imdisp( img_path_nps_all_norm(:,:) );

img_path_all = img_ps2(:,:,1:4:30);
figure; imdisp( img_path_all(:,:), [0.18 0.22] );


%%
close all;

figure
filename = 'nps_100mas_rog.gif';
for n = 1 : 40
    imdisp( img_path_nps_2(:,:,n) );
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
img_path_nps_profile_all = zeros( 128, 40 );
for path_id = 1: 40
    img_path_nps = img_path_nps_2(:,:,path_id);
    for angle = 1:89
        img_path_nps = img_path_nps + imrotate( img_path_nps, angle, 'crop');
    end
    img_path_nps_profile =  img_path_nps(end/2, end/2+1:end ) / 90;
    img_path_nps_profile_all(:,path_id) = img_path_nps_profile / ( max(img_path_nps_profile) );
end

figure;
ix = (0:119) / 120 * ( 1 / 0.7 / 2);
plot( ix, img_path_nps_profile_all(1:120,1:3:20) )
axis([0 0.7 0 1.2])

%%
img_fbp_soft_nps_profile = 0;

img_pwls_nps_profile_2 = 0;
img_pwls_nps_profile_3 = 0;
for angle = 0:365
    img_fbp_soft_nps_profile = img_fbp_soft_nps_profile + imrotate( img_fbp_soft_nps, angle, 'crop');
    img_pwls_nps_profile_2 = img_pwls_nps_profile_2 + imrotate( img_pwls_nps_2, angle, 'crop');
    img_pwls_nps_profile_3 = img_pwls_nps_profile_3 + imrotate( img_pwls_nps_3, angle, 'crop');
end

img_fbp_soft_nps_profile = img_fbp_soft_nps_profile(end/2, end/2+1:end ) / 360;
img_pwls_nps_profile_2 =  img_pwls_nps_profile_2(end/2, end/2+1:end ) / 360;
img_pwls_nps_profile_3 =  img_pwls_nps_profile_3(end/2, end/2+1:end ) / 360;

img_fbp_soft_nps_profile = img_fbp_soft_nps_profile / max( img_fbp_soft_nps_profile );
img_pwls_nps_profile_2 = img_pwls_nps_profile_2 / max( img_pwls_nps_profile_2 );
img_pwls_nps_profile_3 = img_pwls_nps_profile_3 / max( img_pwls_nps_profile_3 );

%%
close all

figure;
plot(ix, img_fbp_soft_nps_profile(1:120), 'r-' , 'LineWidth', 2 );  hold on;
plot(ix, img_pwls_nps_profile_2(1:120), 'g-' , 'LineWidth', 2 );
plot(ix, img_pwls_nps_profile_3(1:120), 'b-' , 'LineWidth', 2 ); 
%plot( ix, img_path_nps_profile_dog(1:120,6:6:40), '--', 'LineWidth', 2 );
plot( ix, img_path_nps_profile_all(1:120,6:5:40), '--', 'LineWidth', 2 );
legend( 'FBP soft kernel', 'PWLS \beta = 2x10^4', 'PWLS \beta = 5x10^4' );
axis([0 0.7 0 1.1])
ylabel('Normalized NPS', 'FontSize', 18 );
xlabel('f_{xy} (mm^{-1})', 'FontSize', 18 );
set( gca , 'FontSize', 14);