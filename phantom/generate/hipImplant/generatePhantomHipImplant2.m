I = zeros( 2400, 1600 ); 
meta.DimSize = size(I) ;
meta.ElementSpacing = [0.16 0.16];

fprintf('Painting tissue patterns... \n');
ellipsesParams = [
       0      0   113   173    0   2;
       0      0   110   170    0   1;
       0     70   70    70     0   2;
       0    -70   70    70     0   2;
       0      0   80    90     0   2;
      60     60   30    70    10   2;
      60    -60   30    70   -10   2;
       
       0    -90   30    16   -60   4;
       0     90   30    16    60   4;
      25    -60   20    12     0   4;
      25     60   20    12     0   4;
     -30    -50   12     8   -50   4;
     -30     50   12     8    50   4;
       0    -90   27    13   -60   3;
       0     90   27    13    60   3;
      25    -60   16     8     0   3;
      25     60   16     8     0   3;
     -30    -50    8     5   -50   3;
     -30     50    8     5    50   3;
     -10    -75   15    15     0   5;
     -10     75   15    15     0   5;
	% line pattern between fillings

	% line pattern away from fillings
	  -2     0   1.1   12     15  1;
	  -6     0   1.1   12     15  1;
	 -10     0   1.1   12     15  1;
	 -14     0   1.1   12     15  1;
	 -18     0   1.1   12     15  1;
	  42     0   1.1   12     15  1;
	  46     0   1.1   12     15  1;
	  50     0   1.1   12     15  1;
	  54     0   1.1   12     15  1;
	  58     0   1.1   12     15  1;];
  
  I = addEllipses(I, meta, ellipsesParams);

writeMetaImage(uint8(I), 'HipImplants.mhd', meta);

figure;
imagesc(I'), colormap gray; axis image;