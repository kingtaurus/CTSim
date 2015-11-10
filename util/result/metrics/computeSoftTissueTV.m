function R = computeSoftTissueTV( img, map )

R = 0;

img = squeeze(img(:,:,ceil(end/2)));

d = abs(img(1:end-1, :) -  img(2:end, :));
R = R + sum( d( map.mapTissue(:)) );

d = abs(img(:, 1:end-1) -  img(:, 2:end));
R = R + sum( d( map.mapTissue(:)) );

R = R / sum( map.mapTissue(:) ) / 2;

fprintf('The total variation of soft tissue is %.2f HU. \n', R );

end