
ellipsesParams = zeros(1, 6 );
ellipsesParams(1,:) = [0    0   100  65  0   1];
ellipsesParamsDensity = [0    0   100  65  0  1];

I = zeros( [512 512], 'single' );
meta2.DimSize = size(I) ;
meta2.ElementSpacing = [0.5 0.5];
I = addEllipses(I, meta2, ellipsesParams);

J = zeros( [512 512], 'single' );
J = addEllipses(J, meta2, ellipsesParamsDensity);

figure;
imagesc(I .* J), colormap gray; axis image;

%%
K = zeros([ 512 512 64], 'single' );
meta.DimSize = size(K) ;
meta.ElementSpacing = [0.5 0.5 1];
density = K;

boneTilteAngle = 0.4 * pi;
boneRaduis = 4.5 / 2;
boneRoationAngle = [ 0.5*pi, 0.25*pi, 0.75*pi, -0.4*pi, 1.4*pi   ];
boneCenterX = [0, -70, 70, -36, 36 ];
boneCenterY = [-50, -28, -28, 45, 45 ];

spineTilteAngle = 0.1 * pi;
spineRaduis = 12.7 / 2;
spineRoationAngle = - 0.2 * pi ;
spineCenterX = 5;
spineCenterY = 38;
                    
airTilteAngle = 0.1 * pi;
airRaduis = 10 / 2;
airRoationAngle = 0.8 * pi ;
airCenterX = 40;
airCenterY = -10; 


for i = 1:meta.DimSize(3)                           
    
    j = round( i - meta.DimSize(3)/2 );
    z = j * meta.ElementSpacing(3);
    
    ellipsesParams = zeros(2, 6);
    ellipsesParams(1,:) = [0    0   100  65   0   1];
                         
    I = zeros( [512 512], 'single' );       
    for t = 1:length(  boneRoationAngle )
        
        ellipsesParams(t +1,:) = [ boneCenterX(t) - z * tan(boneTilteAngle)*sin(boneRoationAngle(t))  ....
            boneCenterY(t) + z * tan(boneTilteAngle)*cos(boneRoationAngle(t))  ...
            boneRaduis ...
            boneRaduis/cos(boneTilteAngle) ...
            boneRoationAngle(t)*180/pi   2];
        
    end          
    
    ellipsesParams(7,:) = [ spineCenterX - z * tan(spineTilteAngle)*sin(spineRoationAngle)  ....
        spineCenterY + z * tan(spineTilteAngle)*cos(spineRoationAngle)  ...
        spineRaduis ...
        spineRaduis/cos(spineTilteAngle) ...
        spineRoationAngle*180/pi   2];
    
      ellipsesParams(8,:) = [ airCenterX - z * tan(airTilteAngle)*sin(airRoationAngle)  ....
        airCenterY + z * tan(airTilteAngle)*cos(airRoationAngle)  ...
        airRaduis ...
        airRaduis/cos(airTilteAngle) ...
        airRoationAngle*180/pi   0];
    
    meta2.DimSize = size(I) ;
    meta2.ElementSpacing = [0.5 0.5];
    I = addEllipses(I, meta2, ellipsesParams);

    imagesc(I .* J), colormap gray; axis image;
 
    pause();
    K(:,:,i) = I;
    density(:,:,i) = J;

end

%%
figure;
imagesc( K(:,:,end/2)); colormap gray; axis image;

%% save to raw data
writeMetaImage(uint8(K), 'HBP1-phantom-3d.mhd', meta);
writeMetaImage(density, 'HBP1-phantom-3d-density.mhd', meta);

