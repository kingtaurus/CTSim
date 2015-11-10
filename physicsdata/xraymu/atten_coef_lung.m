function mus = atten_coef_lung(egy)

mus = zeros(1,5);

compositionData = BodyMaterialCompositionFunc;
mu_lung = body_material_mu( 'lung', egy, compositionData );
mu_muscle = body_material_mu( 'muscle_skeletal', egy, compositionData );
mu_water =  body_material_mu( 'water', egy, compositionData );
mu_skin =  body_material_mu( 'skin', egy, compositionData );

mu_catheter = 2712 *  XrayMu('Al', egy )/10;
mu_lung = 0.4 * mu_lung;

mus(1) = mu_lung;
mus(2) = mu_muscle;
mus(3) = mu_water;
mus(4) = mu_catheter;
mus(5) = mu_skin;

end