function numPhotonsPerPixel = xrayIntensitySimilarDose( photonsPerMm2At1mKeV, noViews, geom )

Pixel2msAt1m = (geom.SDD/1000)^2 / ( geom.detSpacing(1) * geom.detSpacing(end));
numPhotonsPerPixel           = photonsPerMm2At1mKeV * noViews / ( Pixel2msAt1m * geom.noViews ) ;

end