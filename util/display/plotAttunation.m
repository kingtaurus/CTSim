function plotAttunation

e = logspace(1, 4, 100 );

figure;

muSoft = materialAttenuation( e, 'Tissue_Soft_ICRU-44');
muBone = materialAttenuation( e, 'Bone_Cortical_ICRU-44');
muTi = materialAttenuation( e, 'Titanium');
muAm = materialAttenuation( e, 'AK-Amalgam');
muAu = materialAttenuation( e, 'Gold');

loglog( e, muAu, '.-', 'lineWidth', 2 , 'Color', 'k'); hold on;
loglog( e, muAm, ':', 'lineWidth', 2 , 'Color', 'k'); hold on;
loglog( e, muTi, '.', 'lineWidth', 2 , 'Color', 'k'); hold on;
loglog( e, muBone, '-.', 'lineWidth', 2 , 'Color', 'k'); hold on;
loglog( e, muSoft, 'lineWidth', 2 , 'Color', 'k'); hold on;



grid on;
xlabel( 'keV', 'fontSize', 20);
legend( 'Gold', 'Amalgam','Titanium','Bone', 'Soft Tissue'  );

set(gca,'FontSize',20);
axis tight;