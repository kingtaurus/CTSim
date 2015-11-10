
scanAngle = 40;
offsetAngle = 90;

load 'temp.mat';

fileName = 'prostateScans';

phan = loadCTScanPhantom( fileName  );

geom = loadProjectionGeometryCarm( p, scanAngle, offsetAngle );


%% display tomographic angles

tomoAngles = [];
z = -200:1:200;

for  scanAngle = 20 : 10 : 60
    ta = tomoAngleCircluarGeometry( geom, scanAngle, z );
    tomoAngles = [tomoAngles; ta];
end

figure;
plot( z, tomoAngles);
xlabel( 'Vertical distance for rotation axis (mm)' );
ylabel('Full tomograhic angles (degrees)');
legend('20^o','30^o','40^o','50^o','60^o');


%%

scanAngle = 180;
offsetAngle = 90;

p.Geometries.keV.noViews = round( scanAngle / 2 );


geom = loadProjectionGeometryCarm( p, scanAngle, offsetAngle );


sinoAtt = computeAttenuationSinogramFromCTScan( phan, geom );

fprintf('Simulate projection done. \n')

geom.betas = geom.betas;

img = reconFBP( sinoAtt, geom );

fprintf('SAA reconstruction done. \n')

%%

mijshow( img(:,:,130)); 


mijshow( squeeze( img(133, :, :))') ; 
