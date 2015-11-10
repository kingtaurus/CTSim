function imquick( img, windows)

if nargin < 2
    windows = [0.1 0.35];
end

h = figure;
if ndims( img ) == 3
    
    img = rotateSinogram( img, 0, 1 );
    k = round( size(img, 3) / 8 );
    s = ceil( mod( size(img, 3), k ) / 2) + 1;

    imdisp( img(:,:, s : k : end ), windows  );
    imdisp( fliplr( [ squeeze( img(:,round(end*0.55),:))  squeeze( img(:,round(end*0.5),:))  squeeze( img(:,round(end*0.45),:)) ] ), windows );
else
    imdisp( img', windows  );
end

set( h, 'NumberTitle', 'off', 'Name', [inputname(1) ' ' num2str( fix( mod( cputime * 100, 10000) ) )]);
end