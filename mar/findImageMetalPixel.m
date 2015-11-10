function imgMetalMap = findImageMetalPixel( sino, geom, n )

metalThreshold = 1.0;

sino = medianFilterSino( sino, 3 );
imgFBP = reconFBP( sino, geom);
imgMetalMap = imgFBP > metalThreshold;

se =  strel('disk', n);

if length( size( imgMetalMap )) == 2
    imgMetalMap = imdilate( imgMetalMap, se);
else
    for iz = 1:size(imgMetalMap,3)
        imgMetalMap(:,:,iz) =  imdilate( imgMetalMap(:,:,iz), se);
    end
end