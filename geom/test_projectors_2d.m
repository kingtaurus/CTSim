%test_projectors

close all;

n= 128;

I = zeros(n,n);

ix = 1:n;
ix = ix + 0.5;
iix = ones(n,1) * ix;
iiy = iix';

r1 = round(n/12);
r2 = round(n/32);
r3 = round(n/3);

c1x = round(n*2/3);
c1y = round(n/2);

c2x = round(n/2);
c2y = round(n/3);

c3x = round(n/2);
c3y = round(n/2);

I( sqrt((iix-c3x).^2+(iiy-c3y).^2) < r3 ) = 1;
I( sqrt((iix-c1x).^2+(iiy-c1y).^2) < r1 ) = 2;
I( sqrt((iix-c2x).^2+(iiy-c2y).^2) < r2 ) = 5;


img = I;

figure(1)
imagesc(img); colormap gray;


geom.detSpacing = 1;
geom.detSize = 400;
geom.detOffset = 10;


phan.size = size(img);
phan.spacing = [1.5 1.5];
phan.offset  = [10, 10];

geom.reconSize = phan.size;
geom.reconSpacing = phan.spacing;
geom.reconOffset = phan.offset;


geom.SAD = 1000;
geom.ADD = 800;
geom.SDD = geom.SAD + geom.ADD;
geom.noViews = 160;

geom.betas = (0:geom.noViews-1) * 2*pi/geom.noViews;

tic
sino2 = forwardProjectDistanceDrivenMex( img, phan, geom );
toc

tic
sino1 = forwardProjectRayDrivenMex( img, phan, geom );
toc


figure(12)
imagesc(sino1); colormap gray;

figure(13)
imagesc(sino1- sino2); colormap gray;


%%

tic
sinoTF = forwardProjectTrapezoidicFootprintMex( img, phan, geom );
toc



figure(4)
imagesc(sinoTF); colormap gray;


figure(5)
imagesc(sino - sinoTF ); colormap gray;


tic
img = backProjectDistanceDrivenMex( sinoMex, geom );
toc

tic
imgMex = backProjectPixelDrivenMex( sinoMex, geom );
toc

figure(3)
imagesc(imgMex); colormap gray;

figure(4)
imagesc(img); colormap gray;



return;
%%








