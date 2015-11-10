function [muSoftEff, muBoneEff] = getEffectiveAttuation( img, muSoft, muBone )


fsoft = (img <= muSoft * 1.5) & (img > 0.5 * muSoft);
fbone = img > muSoft * 1.5;

muBoneEff = mean(img(fbone(:)));
muSoftEff = mean(img(fsoft(:)));



end