function showPhaserRecons( img, windows, format )

if nargin < 2
    windows = [0.1 0.35];
end

img = rotateSinogram( img, 0, 0 );

figure; imdisp( img(:,:,end:-8:1), windows);
if nargin > 2
    imgname = inputname(1);
    saveas( gcf, [imgname '_all'], format );
end


figure; imdisp( img(:,:,end/2), windows);
if nargin > 2
    imgname = inputname(1);
    saveas( gcf, [imgname '_axail'], format );
end


figure; imdisp( squeeze(img(:,end/2,end:-1:1))' , windows );
if nargin > 2
    imgname = inputname(1);
    saveas( gcf, [imgname '_sagittal'], format );
end


figure; imdisp( squeeze(img(end/2,:,end:-1:1))' , windows );
if nargin > 2
    imgname = inputname(1);
    saveas( gcf, [imgname '_coronal'], format );
end
