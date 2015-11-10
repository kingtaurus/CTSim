function [fsoft, fbone, fmetal, pho] = materialSegmentationMetalCont( img, muSoft, muBone, muMetal )
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
% Meng Wu at Stanford University
% 2012 - 2013


diff = muBone - muSoft;

initailThresholds = [muSoft/2  muSoft+diff/3 muSoft+diff*2/3 muBone+diff/3];

[clusterThresholds, clusterMeans, clusterStds] = segmentationAdaptiveThreshold( img, initailThresholds);

muSoft = clusterMeans(2);
muBone = clusterMeans(4);



stdSoft = clusterStds(2);
stdBone = clusterStds(4);

imgManipulated = img;
imgManipulated( imgManipulated < muSoft + 2 * stdSoft) = muSoft;
imgManipulated( imgManipulated > muBone - 2 * stdBone) = muBone;

fsoft = ( muBone -  imgManipulated) / (muBone - muSoft);

fbone = 1 - fsoft;

fbone( img > 1.2* muBone ) = ( muMetal - img( img > 1.2 * muBone ) ) / (muMetal - 1.2 * muBone);

fbone( fbone < 0 ) = 0;

fmetal = 1 - fbone - fsoft ;

pho = img./ ( muSoft * fsoft + muBone * fbone + muMetal * fmetal );

end