function [fsoft, fbone, fmetal, pho] = materialSegmentationMetalCombinedCont( imgInitial, muSoft, muBone, muMetal )
% Material segmentation with continoust transitision between softtissue and
% bone
%  Input:   
%       imgMAR         - image to segment
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

[~, clusterMeans, clusterStds] = segmentationAdaptiveThreshold( imgInitial, initailThresholds);

muSoft = clusterMeans(2);
muBone = clusterMeans(4);

stdSoft = clusterStds(2);
stdBone = clusterStds(4);

imgManipulated = medfilt3( imgInitial, [3 3 1] );
imgManipulated = imgInitial;
imgManipulated( imgManipulated < muSoft + 2 * stdSoft) = muSoft;
imgManipulated( imgManipulated > muBone - stdBone) = muBone;

fsoft = ( muBone -  imgManipulated) / (muBone - muSoft);
fbone = 1 - fsoft;
fbone( fbone < 0 ) = 0;

fmetal = imgInitial > muMetal ;
fmetal = single( fmetal );
fmetal = imerode( fmetal, [0 1 0; 1 1 1; 0 1 0]);
fmetal = imfilter( fmetal, [0.05 0.1 0.05; 0.1 0.4 0.1; 0.05 0.1 0.05 ]);

fbone =  min( fbone + imfilter( fmetal, [0.05 0.1 0.05; 0.1 0.4 0.1; 0.05 0.1 0.05 ]), 1);


pho =  imgInitial ./ ( muSoft * fsoft + muBone * fbone + muMetal * fmetal );

end