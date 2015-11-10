% function loadGeant4DoseOutput

doseDir = '../CTData/DoseSimData/hipImplant120KV/';
info = dicominfo([doseDir 'center.dcm'] );

noViews = 240;
noFillings = 1;

Geant4ParticleDensity = 1e6;
ReconParticlePerPixel = 8e5;
ReconNoViews          = 480;

pixelSize   = 0.768;
figureName  = 'dose_hip_kv.eps';


Y = dicomread(info);
%figure, imshow(Y, [0 600]);

n = double( info.Rows );
m = double( info.Columns );

doseMap = zeros( [n * m, 1] );

% read output file
for iv = 1:noViews
    
    if mod(iv, 20) == 1
        fprintf( 'Angle %g \n', iv);
    end
    
    for ib = 1 : noFillings
        
        fileName = sprintf( 'dose_%03g_%01g.out', iv, ib );
        
        fid = fopen([doseDir 'output/' fileName ], 'r');
        if fid > 0
            A = fscanf(fid, '%g %g', [2 inf]);
            fclose(fid);
        end
        %index = mod( A(1,:), n * m ) ;
        doseMap( A(1,:)+1 ) = doseMap( A(1,:)+1 ) + A(2,:)';
        
    end
end

%%

SDD = 1500;
HalfDihedralAngles = atan( pixelSize / 2 / SDD );
PixelSolidAngle = 4 * asin( sin(HalfDihedralAngles) * sin(HalfDihedralAngles) );
EquivalentParticlePerPixel  = Geant4ParticleDensity * 180 / pi * PixelSolidAngle;


ScaledDoseMap = doseMap * ReconParticlePerPixel * ReconNoViews / ( EquivalentParticlePerPixel * noViews);

ScaledDoseMap = reshape( ScaledDoseMap, [m n] );
doseMapFilt = medfilt2( ScaledDoseMap, [5 5] );

%% 

map.mapTissue = imerode( ( Y' == 100 ) | ( Y' == 200 ), ones(5));
map.mapBone = imerode( ( Y' == 300 ) | ( Y' == 400 ), ones(5) );
map.mapFillings = ( Y' == 500 ) | ( Y' == 600 );
computeMaterialDoes( map, doseMapFilt );
doseMapFilt( Y' == 0 ) = 0;
 
 
%% display

[phan, map, roi] = loadMaterialsPhantom(p);

y = ((-(phan.size(1)-1)/2: (phan.size(1)-1)/2) + phan.offset(1)) * phan.spacing(1);
x = ((-(phan.size(2)-1)/2: (phan.size(2)-1)/2) + phan.offset(2)) * phan.spacing(2);

[xx, yy] = meshgrid(x, y);
iy = ((-(phan.reconSize(1)-1)/2: (phan.reconSize(1)-1)/2) + phan.reconOffset(1)) * phan.reconSpacing(1);
ix = ((-(phan.reconSize(2)-1)/2: (phan.reconSize(2)-1)/2) + phan.reconOffset(2)) * phan.reconSpacing(2);

[ixx, iyy] = meshgrid(ix, iy);

imgDoseMap = interp2( xx, yy, doseMapFilt, ixx, iyy );

figure; 
imagesc( flipud( imgDoseMap' ), [0 0.1]  );
colormap( hot(256));
axis xy; axis off;
h = colorbar;
title(h,'Gy','FontSize',16)
set(gca,'FontSize',16);
saveas(gcf, figureName,'jpg');
