function [fsoft, fbone, pho] = materialSegmentationCont( imgInitial, muSoft, muBone )
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


diff = muBone - muSoft;

initailThresholds = [muSoft/2  muSoft+diff/3 muSoft+diff*2/3 muBone+diff/3];

[~, clusterMeans, clusterStds] = segmentationAdaptiveThreshold( imgInitial, initailThresholds);

muSoft = clusterMeans(2);
muBone = clusterMeans(4);

stdSoft = clusterStds(2);
stdBone = clusterStds(4);

imgManipulated = medfilt3( imgInitial, [3 3 1] );
imgManipulated( imgManipulated < muSoft + 2 * stdSoft) = muSoft;
imgManipulated( imgManipulated > muBone - 2 * stdBone) = muBone;

fsoft = ( muBone -  imgManipulated) / (muBone - muSoft);

fbone = 1 - fsoft;

pho = imgInitial./ ( muSoft * fsoft + muBone * fbone) ;


end
