function phan = loadCTScanPhantom( fileName  )


[phan.phantomMaterials, metaMaterials] = readMetaImage( [ fileName '.mhd']);

phan.size           = metaMaterials.DimSize; % number of pixels in ground truth materials image
phan.spacing        = metaMaterials.ElementSpacing; % input spacing (mm)
phan.offset         = zeros( size(phan.size) );

phan.reconSize      = round(phan.size);
phan.reconSpacing   = phan.spacing ;
phan.reconOffset    = phan.offset ;

end