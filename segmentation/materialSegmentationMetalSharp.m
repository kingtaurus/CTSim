function [fsoft, fbone, fmetal, pho] = materialSegmentationMetalSharp( u, muSoft, muBone, muMetal )
% Material segmentation with sharp transitision between softtissue and
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

r1 = 1.3;
fsoft = u <= muSoft * r1;
fmetal =  u >= muBone * r1;
fbone = ~fsoft & ~fmetal ;

pho = zeros( size(u) );
pho( fsoft )    = u( fsoft ) / muSoft;
pho( fbone )    = u( fbone ) / muBone;
pho( fmetal )   = u( fmetal ) / muMetal;

end