close all;

load 'burn01.mat';

[optimizer, metric] = imregconfig('monomodal');

imgAirReg = imgAir;

for i = 1 : size( imgAir, 3 )
    
    slice = imregister(imgAir(:, :, i), imgKr(:,:,i), 'rigid', optimizer, metric);
    
    imgAirReg( :,:,i) = slice;
    
end



%%
solid = imgAirReg > 0.2;

final = solid;

solid_blurred = solid; 

for i = 1 : size( imgAir, 1 )
    solid_blurred( i, :, : ) = imdilate( solid(i,:,:) , ones(2) ); 
end

final = final | solid_blurred;

for i = 1 : size( imgAir, 1 )
    solid_blurred( :, i, : ) = imdilate( solid(:,i,:) , ones(2) ); 
end

final = final | solid_blurred;

for i = 1 : size( imgAir, 3 )
    solid_blurred( :, :, i ) = imdilate( solid(:,:,i) , ones(2) ); 
end


final = final | solid_blurred;

final( :, :, 395 : end ) = false;


%% 

imgSub = imgKr - imgAirReg;

%%
x = [95 182];
y = [61 178];

att_curve = zeros( 1, size( imgAir, 3));

for i = 1 : length( att_curve )
    
    slice = imgSub(x(1):x(2), y(1):y(2),i);
    valid =  ~final( x(1):x(2), y(1):y(2),i);
    
    att_curve(i) = mean( slice( valid(:) ) );
    
end

figure; plot( att_curve ); 

figure;
imagesc( squeeze( imgSub(:,end/2,:) )' , [-0.01 0.02] ); 


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

