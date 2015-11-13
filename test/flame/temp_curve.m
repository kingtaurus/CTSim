close all;

% of course, we need subtraction images
imgSub = imgKr - imgAir;

%% Get image pixel that are not gas

% segmentation
solid = imgAir > 0.18;
solid = solid | imgKrReg > 0.18 ;

% mophological blurring in all dimension
final = solid;
solid_blurred = solid; 

% in x direction
for i = 1 : size( imgAir, 1 )
    solid_blurred( i, :, : ) = imdilate( solid(i,:,:) , ones(3) ); 
end
final = final | solid_blurred;

% in y direction
for i = 1 : size( imgAir, 1 )
    solid_blurred( :, i, : ) = imdilate( solid(:,i,:) , ones(3) ); 
end
final = final | solid_blurred;

% in z direction
for i = 1 : size( imgAir, 3 )
    solid_blurred( :, :, i ) = imdilate( solid(:,:,i) , ones(3) ); 
end
final = final | solid_blurred;

final = final | abs( imgSub ) > 0.04;

% for slices that have porous media below system resolution
% final( :, :, 241 : end ) = false;


%% Now let's compute average density for each slice
close all;

% bounding box with the
x = [77 180];
y = [75 180];

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

