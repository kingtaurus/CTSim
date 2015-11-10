function [fsoft, fbone, pho] = materialSegmentationFree( img, muSoft, muBone )
% Material segmentation with continoust transitision between softtissue and
% bone
%  Input:   
%       img         - image to segment
%       muSoft      - attenuation of soft tissue
%       muBone      - attenuation of bone
%  Output:
%       fsoft       - segmentation of soft tissue
%       fbone       - segmentation of bone
%       pho         - denisty
%
% Meng Wu at Stanford Universityiz
% 2012 - 2013

img( img < muSoft ) = muSoft;

fsoft = ( muSoft ./ img ) ;
fsoft( fsoft >= 1 ) = 1;
fbone = 1 - fsoft;

pho = img./ ( muSoft * fsoft + muBone * fbone) ;


end
