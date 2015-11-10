function showBinnigErrorMap( spectrum, BinSizes, Phis, Thetas, w )
% show error map of spectrum binning for different thickness of materials


n = 64;
BinSizes    = BinSizes(:);
Phis        = Phis(:);
Thetas      = Thetas(:);

photonsPerEnergyBin     = spectrum.photonsPerEnergyBinOriginal * spectrum.DQE;
[phis, thetas, softThickness, boneThickness] = generateMaterialThicknessMap( n, spectrum );
[ phiSoft, thetaSoft, muEffective,  PhiAll, ThetaAll ] = attenuationCoefficientDecompose( 'Tissue_Soft_ICRU-44', spectrum );

oc = ones( size(thetas) );

yTruth =  log(  sum ( (photonsPerEnergyBin  * oc') .*  exp( - PhiAll * phis' - ThetaAll * thetas' ) ) )';

yApprox =  log(  sum ( (BinSizes  * oc') .*  exp( - Phis * phis' - Thetas * thetas' ) ) )';

errorMap = reshape( yTruth - yApprox, [n, n] );


figure;
imagesc( softThickness, boneThickness, errorMap, w ); axis xy;
colorbar; 
ylabel ('Bone Thickness (cm)', 'fontsize', 22);
xlabel ('Soft Tissue Thickness (cm)', 'fontsize', 22);
set( gca, 'fontsize', 22);
