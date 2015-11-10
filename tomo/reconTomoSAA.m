function img = reconTomoSAA( sino, geom )


img = backProjectPixelDrivenMex(sino, geom );

w =  backProjectPixelDrivenMex( ones( size(sino), 'single'), geom );

img = img./w;


end