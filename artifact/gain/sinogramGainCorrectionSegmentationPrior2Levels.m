function [sino] = sinogramGainCorrectionSegmentationPrior2Levels( sino, geom, imgPrior, secondSinoIndex, level, softThreshold )
% [sino] = sinogramGainCorrectionSegmentationPrior( sino, geom, redundantPixel1, redundantPixel2, imgPrior )
% note:
%   degree - degree of polynomial fit
%
%
% Meng Wu @ Stanford Univercity
% 2014.2

if nargin < 6
    softThreshold = false;
end
fprintf( 'Sinogram gain linaer correction using redundant data : \n');

nu = geom.detSize(1);
du = geom.detSpacing(1);
ou = geom.detOffset(1);

noViews = geom.noViews;
SDD = geom.SDD;
SAD = geom.SAD;
betas = geom.betas - min(geom.betas);

% array of detector positions in mm
u = ( ( -(nu-1)/2:(nu-1)/2) + ou ) * du ;
gamma = atan( - u / SDD );
Tau  = abs(geom.betas(end) - geom.betas(1)) - pi;

weightMap = zeros( geom.detSize(1), noViews );
beta = zeros( geom.detSize(1), noViews );
theta = zeros( geom.detSize(1), noViews );

% label pixels in the fan beam geomerty
for iv = 1:noViews
    weightMap(:,iv) =  shortScanSliverWeight( gamma, betas(iv),  Tau);
    beta(:,iv) =  betas(iv);
    theta(:,iv) = + gamma ;
end

% forward projection segemenated proir image to get a estimate of total
% thicknessProjections
%level = graythresh(imgPrior);
fprintf('\tSegmentation level = %2.3f, %2.3f. \n', level(1), level(2) );


if softThreshold
    imgBW1 = imgPrior ./ level(1);
    imgBW1( imgPrior < level(1) ) = 0 ;
    
    imgBW2 = imgPrior ./ level(2);
    imgBW2( imgPrior < level(2) ) = 0 ;
    
else
    imgBW1 = imgPrior > level(1);
    imgBW2 = imgPrior > level(2);
end
thicknessTissue = forwardProjectMex( single(imgBW1),geom );
thicknessBone = forwardProjectMex( single(imgBW2),geom );

% figure; imdisp( squeeze(thicknessTissue(180,1:4:end,:) ) )
% figure; imdisp( squeeze(thicknessBone(180,1:4:end,:) ) )

% find coresponding redundant pixels and total thicknessProjections
redundantSino1 = weightMap < 0.9 & weightMap > 0.1;
redundantSino1(:,secondSinoIndex:end) = 0;
% figure;
% imagesc( weightMap ); xlabel 'beta', ylabel 'gamma'
beta2 = beta(redundantSino1) + pi - 2 * theta(redundantSino1 );
gamma2 = - theta(redundantSino1 );
[bb, gg] = meshgrid( betas(secondSinoIndex:end), gamma );

k = 5;
if ndims( sino ) == 3
    
    numberOfRedundentPixels = sum(redundantSino1(:));
    redundantPixel1 = zeros( numberOfRedundentPixels * k, 1);
    redundantPixel2 = zeros( numberOfRedundentPixels * k, 1);
    redundantThicknessTissue = zeros( numberOfRedundentPixels * k, 1);
    redundantThicknessBone = zeros( numberOfRedundentPixels * k, 1);
    
    for i = 0:k-1
        centralRow = squeeze( sino(end/2 + i - (k-1)/2,:,:));
        centralRowThickness = squeeze( thicknessTissue(end/2 + i - (k-1)/2,:,:));
        redundantThicknessTissue(i*numberOfRedundentPixels+1 : (i+1)*numberOfRedundentPixels) = centralRowThickness( redundantSino1(:) );
        centralRowThickness = squeeze( thicknessBone(end/2 + i - (k-1)/2,:,:));
        redundantThicknessBone(i*numberOfRedundentPixels+1 : (i+1)*numberOfRedundentPixels) = centralRowThickness( redundantSino1(:) );
        
        
        redundantPixel1(i*numberOfRedundentPixels+1 : (i+1)*numberOfRedundentPixels) = centralRow( redundantSino1(:) );
        redundantPixel2(i*numberOfRedundentPixels+1 : (i+1)*numberOfRedundentPixels) = interp2( bb, gg, centralRow(:,secondSinoIndex:end), beta2, gamma2 );
        
    end
else
    redundantPixel1 = sino( redundantSino1(:) );
    redundantThicknessTissue = thicknessTissue( redundantSino1(:) );
    redundantThicknessBone = thicknessBone( redundantSino1(:) );
    redundantPixel2 = interp2( bb, gg, sino(:,secondSinoIndex:end), beta2, gamma2 );
    
end

redundantPixel1 = redundantPixel1( ~isnan( redundantPixel2 ));
redundantThicknessTissue = redundantThicknessTissue(~isnan( redundantPixel2 ));
redundantThicknessBone = redundantThicknessBone(~isnan( redundantPixel2 ));
redundantPixel2 = redundantPixel2( ~isnan( redundantPixel2 ));

fprintf( '\tNumber of redundant pixels: %i \n', length(redundantPixel2));
fprintf( '\tPercentage of total pixels: %2.3f \n', length(redundantPixel2)/numel(sino) );

% linear regression
n = length(redundantPixel2);
X = [ redundantPixel2 redundantThicknessTissue redundantThicknessBone ones(n, 1) ];
p = ( X' * X )^(-1) * X' * redundantPixel1;
Z = ( X' * X )^(-1);

fprintf('\tModel 2: sino1 = %2.2f * sino2 +  %2.2f * thickness1 +  %2.2f * thickness2 +  %2.2f \n', p(1), p(2), p(3) , p(4) );
residues = redundantPixel1 - ( p(1) * redundantPixel2 + p(2) * redundantThicknessTissue  + p(3) * redundantThicknessBone + p(4) );
sigma =  sqrt( sum( residues.^2 ) / (n-1));
for i = 1:4
    z = p(i) / ( sigma * sqrt( Z(i,i) ) ) ;
    fprintf('\tp(%i) = %2.3f, \tt-test = %2.3f, \tp-value = %2.3e. \n', i,  p(i), z, 1-tcdf(abs(z), n-4) );
end

fprintf('\tStandard Error = %2.3f, R-squared = %2.4f\n', sigma, 1 - sigma^2 / std(redundantPixel1) );


% figure;
% subplot(221);
% plot(redundantPixel2, redundantPixel1, '.' ); xlabel 'mv', ylabel 'kv';
% subplot(222);
% plot(redundantThicknessTissue, redundantPixel1 - redundantPixel2 * p(1) - p(4), '.' ); xlabel 'tissue', ylabel 'sub';
% subplot(223);
% plot(redundantThicknessBone, redundantPixel1 - redundantPixel2 * p(1) - redundantThicknessTissue * p(3) - p(4), '.' ); xlabel 'bone', ylabel 'sub';
% subplot(224);
% plot(redundantPixel1, residues, '.r' ); xlabel 'kv', ylabel 'residue'; axis tight;

if ndims( sino ) == 3
    sino(:,:,secondSinoIndex:end) = p(1) * sino(:,:,secondSinoIndex:end) ...
        + p(2) * thicknessTissue(:,:,secondSinoIndex:end)  + p(3) * thicknessBone(:,:,secondSinoIndex:end) + p(4);
else
    sino(:,secondSinoIndex:end) = p(1) * sino(:,secondSinoIndex:end) ...
        + p(2) * thicknessTissue(:,secondSinoIndex:end) + p(3) * thicknessBone(:,secondSinoIndex:end) + p(4);
end

fprintf('Done!\n\n\n');
