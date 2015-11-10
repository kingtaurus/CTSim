function [ imgGtAtt, imgGtHu ] = computeGroundTruth(phan, spectrum)
% Compute ground truth image both in 2D and 3D (only central slice for 3D)
% input:
%       phan        - phantom parameters
%       spectrum    - spectrum parameters
% output:
%       imgGtAtt - ground truth images in Att
%       imgGtHu -  ground truth images in Hu
%
% Meng Wu
% 2013.1

y = ((-(phan.size(1)-1)/2: (phan.size(1)-1)/2) + phan.offset(1)) * phan.spacing(1);
x = ((-(phan.size(2)-1)/2: (phan.size(2)-1)/2) + phan.offset(2)) * phan.spacing(2);

[xx, yy] = meshgrid(x, y);

iy = ((-(phan.reconSize(1)-1)/2: (phan.reconSize(1)-1)/2) + phan.reconOffset(1)) * phan.reconSpacing(1);
ix = ((-(phan.reconSize(2)-1)/2: (phan.reconSize(2)-1)/2) + phan.reconOffset(2)) * phan.reconSpacing(2);

[ixx, iyy] = meshgrid(ix, iy);


%% Monochromatic image
imgGtPhanSizeAtt = materialImageToAttenuation(phan, spectrum);

%imgGtPhanSizeAtt = single( phan.phantomMaterials );

if length(size(imgGtPhanSizeAtt)) == 3
    noSlices = size(imgGtPhanSizeAtt, 3);
    
    if phan.reconSpacing(3) < 2 * phan.spacing(3)
        
        imgGtAtt = interp2( xx, yy, squeeze(imgGtPhanSizeAtt( :,:,ceil(noSlices/2))), ixx, iyy );
        
    else
        imgGtAtt = interp2( xx, yy, sum(imgGtPhanSizeAtt( :,:,ceil(noSlices/2)-1:ceil(noSlices/2)+2),3)/4 , ixx, iyy );
        
    end
else
    imgGtAtt = interp2( xx, yy, imgGtPhanSizeAtt, ixx, iyy );
end

imgGtAtt( isnan(imgGtAtt) ) = 0;

% convert attenuation coefficients to Hounsfield Unit
imgGtHu = convertMonoAttToHu( imgGtAtt, spectrum);


end