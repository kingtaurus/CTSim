function tv = computeTV( img, map )

R = 0;

d = abs(img(1:end-1, :) -  img(2:end, :));
R = R + sum( d( map(:)) );

d = abs(img(:, 1:end-1) -  img(:, 2:end));
R = R + sum( d( map(:)) );


tv  = R / sum( map(:) ) / 2;

end