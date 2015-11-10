function [ fsoft, fbone1, fbone2, fAl, fCa ] =  materialSegmentationMultBHC( u0, muBone1, muBone2, muAl, muCa )


u0 = medfilt2( u0, [5 5]);

fCa = single( u0 > 0.86 );
%fCa = imdilate( fCa, ones(3));
fAl = single( u0 > 0.75 & u0 <= 0.86 );
fbone2 = single( u0 > 0.60 & u0 <= 0.75 );
fbone1 = single( u0 > 0.40 & u0 <= 0.60 );
fsoft = 1 - (fCa + fAl + fbone1 + fbone2 );



figure;
imagesc( fsoft + 2 * fbone1 + 3 * fbone2 + 4 * fAl + 5 * fCa );


end