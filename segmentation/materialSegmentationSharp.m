function [fsoft, fbone, pho] = materialSegmentationSharp( img, muSoft, muBone )
% Material segmentation between softtissue and
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
% Meng Wu at Stanford University
% 2012 - 2013

r1 = 1.5;

fsoft = img <= muSoft * r1;
fbone = img > muSoft * r1;

pho = zeros( size( img ) );
pho( fsoft ) = img( fsoft ) / muSoft;
pho( fbone ) = img( fbone ) / muBone;



end