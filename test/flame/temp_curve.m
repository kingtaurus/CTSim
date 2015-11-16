clear all;
close all;

data = load( 'burn_01.mat');
imgKr = data.img;

data = load( 'air_hot_02.mat');
imgAir = data.img;

clear data;

%% image registration, you may skip this step

imgKrReg = imgKr;

[optimizer, metric] = imregconfig('monomodal');
for i = 1 : size( imgAir, 3 )
    if mod(i, 50) == 0
        fprintf('(%i/%i)... ', i, size(imgAir, 3 ) );
    end
    slice = imregister(imgKr(:,:,i), imgAir(:, :, i), 'rigid', optimizer, metric);
    imgKrReg( :,:,i) = slice;
end

figure(23); imdisp( imgKr(:,end/2,:) - imgAir(:,end/2,:)  , [-0.1 0.1] );

%% Get image pixel that are not gas

imgSub = imgKr - imgAir;

% segmentation
solid = imgAir > 0.22;
solid = solid | imgKr > 0.22;

% mophological blurring in all dimension
final = solid;
solid_blurred = solid; 

% in x direction
for i = 1 : size( imgAir, 1 )
    solid_blurred( i, :, : ) = imdilate( solid(i,:,:) , ones(5) ); 
end
final = final | solid_blurred;

% in y direction
for i = 1 : size( imgAir, 1 )
    solid_blurred( :, i, : ) = imdilate( solid(:,i,:) , ones(5) ); 
end
final = final | solid_blurred;

% in z direction
for i = 1 : size( imgAir, 3 )
    solid_blurred( :, :, i ) = imdilate( solid(:,:,i) , ones(5) ); 
end
final = final | solid_blurred;

final = final | abs( imgSub ) > 0.04;

% for slices that have porous media below system resolution
 final( :, :, 370 : end ) = false;


%% Now let's compute average density for each slice
close all;

% bounding box with the
x = [123 284];
y = [180 290];

%final( 120:140,120:140,:) = true;

att_curve = zeros( 1, size( imgAir, 3));

for i = 1 : length( att_curve )
    
    slice = imgSub(x(1):x(2), y(1):y(2),i);
    valid =  ~final( x(1):x(2), y(1):y(2),i);
    
    att_curve(i) = mean( slice( valid(:) ) );
    
end

figure; plot( att_curve ); 
xlabel 'slice #', ylabel 'attenuation';

figure;
imagesc( squeeze( imgSub(:,end/2,:) )' , [-0.01 0.02] ); axis image

figure;
slice = imgSub(:,end/2,:);
slice( final(:,end/2,:) ) = 0;
imagesc( squeeze( slice )' , [-0.01 0.02] ); axis image


%%
figure;
slice = imgSub(:,:,400);
slice( final(:,:,400) ) = 0;

imagesc(slice , [-0.01 0.02] ); 

figure;
slice = imgSub(:,:,375);
slice( final(:,:,375) ) = 0;

imagesc(slice , [-0.01 0.02] ); 

figure;
slice = imgSub(:,:,350);
slice( final(:,:,350) ) = 0;

imagesc(slice , [-0.01 0.02] );

