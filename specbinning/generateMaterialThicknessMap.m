function [phis, thetas, thickness1, thickness2] = generateMaterialThicknessMap( n, spectrum, ...
    material1, material2,  maxThickness1, maxThickness2 )
% samples material thickness uniformly over the ranges, and the

if nargin < 3
    material1 = 'Tissue_Soft_ICRU-44';
    material2 = 'Bone_Cortical_ICRU-44';
end

if nargin < 5
    maxThickness1 = 40;
    maxThickness2 = 10;
end

[ phi1, theta1 ] = attenuationCoefficientDecompose( material1, spectrum );
[ phi2, theta2 ] = attenuationCoefficientDecompose( material2, spectrum );

thickness1 = linspace(0, maxThickness1, n );
thickness2 = linspace(0, maxThickness2, n );

[ thicknessMap1, thicknessMap2 ] = meshgrid( thickness1, thickness2 );

% create the array of thetas and phis
phis = thicknessMap1 * phi1 + thicknessMap2  * phi2;
thetas = thicknessMap1 * theta1 + thicknessMap2  * theta2;

phis = phis(:);
thetas = thetas(:);