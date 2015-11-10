function  [ sinoOut, geomOut ] = truncationCorrectionTwoPass(  sinoIn, geomIn, mu, edgeWidth, extendWidth, reconFunc )
% function  [ sinoOut, geomOut ] = truncationCorrectionTwoPass(  sinoIn, geomIn, mu, edgeWidth, extendWidth, reconFunc )
% Two-passes truncation correction alogrithm for CT
% 
% Based on SIEMENS white paper:
% Gysbrechts et al., "HD Field of View in Computed Tomography for RT Planning"
%
% Meng Wu @ Stanford University
% 2014

fprintf('Two-pass truncation artifact correction and extended FOV algorithm: \n');

fprintf('\t1) first pass truncation correction ... \n');
nu = geomIn.detSize(1);
[ sinoOut, geomOut, truncatedView ]= truncationCorrectionWaterCylinderFittingTotalConsistency( sinoIn, geomIn, mu, edgeWidth, extendWidth );

fprintf('\t2) image reconstruction ... \n');
img = reconFunc( sinoOut, geomOut );


fprintf('\t3) generate synthetic image ... \n');
x = ( -(geomIn.reconSize(1)-1)/2:(geomIn.reconSize(1)-1)/2);
y = ( -(geomIn.reconSize(2)-1)/2:(geomIn.reconSize(2)-1)/2);
[xx, yy] = meshgrid( y, x);
outSideMap = sqrt( xx.^2 + yy.^2 ) > ( geomIn.FOV / geomIn.reconSpacing(1) / 2 - 4 );

for iz = 1 : size( img, 3)    
    slice = img(:, :, iz );
    slice( outSideMap & slice < mu / 2 ) = 0;
    slice( outSideMap & slice > mu / 2 & slice < mu * 1.5 ) = mu ;
    img(:, :, iz ) = slice ;
    
end


fprintf('\t4) compute synthetic projections ... \n');
% compuate the sinogram error
sinoSynthetic = forwardProjectMex( img, geomOut );

%sinoquick( sinoSynthetic );

fprintf('\t5) second pass truncation correction ... \n');

sinoOut = truncationCorrectionWaterCylinderFitting( sinoIn - sinoSynthetic(:,extendWidth+1: extendWidth+nu, :), geomIn, mu, edgeWidth, extendWidth, truncatedView );

%sinoquick( sinoOut );

fprintf('\t6) combine synthetic projection with truncation correction ... \n');

sinoOut = sinoOut + sinoSynthetic;
sinoOut(:,extendWidth+1: extendWidth+nu, :) = sinoIn;

fprintf('2-Pass Truncation Correction Done. \n');

end