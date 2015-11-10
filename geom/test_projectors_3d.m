%test_projectors


n= 256;

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

%img(c2x,c2y) = 1;

figure(1)
imagesc(I); colormap gray;

img = zeros( [n, n, n/2], 'single');

for i = 1:n/2
   img(:,:,i) = I ;
end
img(:,:,n / 8 * 3 :end) =  img(:,:,n/8*3 :end) / 2 ;


geom.detSpacing = [0.7 0.58];
geom.detSize = [512 256];
geom.detOffset = [ 10 0];

phan.size = [n, n, n/2];
phan.spacing = [0.5, 0.5, 0.5] * 2;
phan.offset  = [10, 10, 0];

geom.reconSize = phan.size;
geom.reconSpacing = phan.spacing;
geom.reconOffset = phan.offset;

geom.SAD = 1000;
geom.ADD = 500;
geom.SDD = geom.SAD + geom.ADD;
geom.noViews = 120;

geom.betas = (0:geom.noViews-1) * 2*pi/geom.noViews;
geom.couchZ = zeros(1, geom.noViews);

geom.flatPanel = 0;
geom.map = true( [n n] );

%%

tic
sino2 = cbct_sf_cuda( 'proj,tf', int32(geom.reconSize), single(geom.reconSpacing), ...
        single(geom.reconOffset), int32(geom.detSize), single(geom.detSpacing), ...
        single(geom.detOffset), single(geom.SAD), single(geom.ADD), single(geom.SDD), ...
        int32(geom.noViews), double(geom.betas), single(geom.couchZ), int32(geom.flatPanel), single(img), logical(geom.map) );
toc

figure(1)
imagesc(squeeze(  sino2(end/2,:,:) ) ); colormap gray;

figure(2)
imagesc(squeeze(  sino2(:,:,end/2) ) ); colormap gray;

figure(3)
imagesc(squeeze(  sino2(:,end/2,:) ) ); colormap gray;


%%

tic
sino = cbct_geom_mex( 'proj,dd', int32(geom.reconSize), single(geom.reconSpacing), ...
        single(geom.reconOffset), int32(geom.detSize), single(geom.detSpacing), ...
        single(geom.detOffset), single(geom.SAD), single(geom.ADD), single(geom.SDD), ...
        int32(geom.noViews), double(geom.betas), single(geom.couchZ), int32(geom.flatPanel), single(img), logical(geom.map) );
toc

tic
sino1 = cbct_geom_mex( 'proj,tf', int32(geom.reconSize), single(geom.reconSpacing), ...
        single(geom.reconOffset), int32(geom.detSize), single(geom.detSpacing), ...
        single(geom.detOffset), single(geom.SAD), single(geom.ADD), single(geom.SDD), ...
        int32(geom.noViews), double(geom.betas), single(geom.couchZ), int32(geom.flatPanel), single(img), logical(geom.map) );
toc


tic
sino2 = cbct_sf_cuda( 'proj,tf', int32(geom.reconSize), single(geom.reconSpacing), ...
        single(geom.reconOffset), int32(geom.detSize), single(geom.detSpacing), ...
        single(geom.detOffset), single(geom.SAD), single(geom.ADD), single(geom.SDD), ...
        int32(geom.noViews), double(geom.betas), single(geom.couchZ), int32(geom.flatPanel), single(img), logical(geom.map) );
toc

%%

close all
 figure(1)
 imagesc(squeeze( sino1(:,:,end/2) ) ); colormap gray;

figure(2)
imagesc(squeeze(  sino2(:,:,end/2) ) ); colormap gray;

figure(3)
imagesc(squeeze( sino2(end/2,:,:) - sino1(end/2,:,:) ) ); colormap gray;

figure(4)
imagesc(squeeze( sino2(:,:,10) - sino1(:,:,10) ) ); colormap gray;


max( abs( sino2(:) - sino1(:)))
%%


tic
img1 = cbct_geom_mex( 'back,tf', int32(geom.reconSize), single(geom.reconSpacing), ...
        single(geom.reconOffset), int32(geom.detSize), single(geom.detSpacing), ...
        single(geom.detOffset), single(geom.SAD), single(geom.ADD), single(geom.SDD), ...
        int32(geom.noViews), double(geom.betas), single(geom.couchZ), int32(geom.flatPanel), single(sino2), logical(geom.map) );
toc

figure(5)
imagesc(squeeze(img1(:,:,end/2))); colormap gray;


%%


tic
img2 = cbct_sf_cuda( 'back,tf', int32(geom.reconSize), single(geom.reconSpacing), ...
        single(geom.reconOffset), int32(geom.detSize), single(geom.detSpacing), ...
        single(geom.detOffset), single(geom.SAD), single(geom.ADD), single(geom.SDD), ...
        int32(geom.noViews), double(geom.betas), single(geom.couchZ), int32(geom.flatPanel), single(sino2), logical(geom.map) );
toc

figure(6)
imagesc(squeeze(img2(:,:,end/2))); colormap gray;


%%

figure(7)
imagesc(squeeze(img2(:,:,end/2)- img1(:,:,end/2) )); colormap gray;


figure(8)
imagesc(squeeze(img2(:,:,1)- img1(:,:,1) )); colormap gray;


