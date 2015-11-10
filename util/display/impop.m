function impop( img, windows, format )


if nargin < 2
    windows = [0.1 0.35];
end

h = figure;
if ndims( img ) == 3
    imdisp( img(:,:,end/2), windows  );
else
    imdisp( img, windows  );
end
set( h, 'NumberTitle', 'off', 'Name', inputname(1));

if nargin > 2
    imgname = inputname(1);
    saveas( gcf, imgname, format );
end



end