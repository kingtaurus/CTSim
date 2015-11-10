function sigma = computeSoftTissueStd( img, map )

img = squeeze(img(:,:,ceil(end/2)));

sigma = std( img( map.mapTissue(:) ) );

fprintf('The standard deviation of soft tissue is %.2f HU. \n', sigma );

end