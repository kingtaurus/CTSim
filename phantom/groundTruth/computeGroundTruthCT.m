function [ imgGtAtt, imgGtHu ] = computeGroundTruthCT(phan, spectrum)


y = [ -(phan.size(1)-1)/2: (phan.size(1)-1)/2] * phan.spacing(1) + phan.offset(1);
x = [ -(phan.size(2)-1)/2: (phan.size(2)-1)/2] * phan.spacing(2) + phan.offset(2);

[xx, yy] = meshgrid(x, y);

iy = [ -(phan.reconSize(1)-1)/2: (phan.reconSize(1)-1)/2] * phan.reconSpacing(1) + phan.reconOffset(1);
ix = [ -(phan.reconSize(2)-1)/2: (phan.reconSize(2)-1)/2] * phan.reconSpacing(2) + phan.reconOffset(2);

[ixx, iyy] = meshgrid(ix, iy);



%covert material phantom
imGtAttPhanSize = size( phan.phantomMaterials );

for m = 1 : length( phan.materialIndexes )

    [~, materialAtt] = materialAttenuation( spectrum.energyAverage, phan.materialNames{m}, spectrum.photonsTotal );
    
    imGtAttPhanSize( phan.phantomMaterials == phan.materialIndexes{m} ) = materialAtt;
    
end

imGtAttPhanSize = reshape( imGtAttPhanSize, phan.size);

if ndims( imGtAttPhanSize ) == 3
imGtAttPhanSize = imGtAttPhanSize( :,:, round( end/2 ) );
end

imgGtAtt = interp2( xx, yy, imGtAttPhanSize, ixx, iyy );
imgGtAtt( isnan(imgGtAtt) ) = 0;

% convert attenuation coefficients to Hounsfield Unit
imgGtHu = convertMonoAttToHu( imgGtAtt, spectrum);

end