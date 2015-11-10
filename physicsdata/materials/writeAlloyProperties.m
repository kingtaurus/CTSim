%% Get Components' Properties

% define materials and read NIST values
materialNames = {'Mercury', 'Silver', 'Tin', 'Copper'};
noMaterials = length(materialNames);
readParameters; % for getting materialsDir
[E_comp, mu_comp, rho_comp, mu_by_rho_comp, mu_by_rho_en_comp] = readMaterialAttenuations(materialNames, materialsDir);

% define composition by mass (reasonable guesstimates for these ratios for Amalgam were made by reading through the Wikipedia article on Amalgam)
m_comp = [0.50, 0.27, 0.13, 0.10]; % in g


%% Compute Components' Properties

% densities
rho_comp = [rho_comp{:}]; % in g/cm^3

% volumes
V_comp = m_comp ./ rho_comp;


%% Alloy Values

% name for the alloy
alloyName = 'AK-Amalgam';

% mass is sum of masses
m_alloy = sum(m_comp);

% volume is sum of volumes (approximately, does not take into account for volume change effects due to packing in an alloy)
V_alloy = sum(V_comp);

% density is mass/volume
rho_alloy = m_alloy / V_alloy;

% energy bins at which the attenuation values are computed
E_alloy =  unique(sort([E_comp{:}]));

% pre-allocate arrays for attenuation values at all energy bins
mu_alloy = zeros(size(E_alloy));
mu_by_rho_alloy = zeros(size(E_alloy));
mu_by_rho_en_alloy = zeros(size(E_alloy));

% compute attenuation values at all energy bins
for e = 1:length(E_alloy)
	% get all components' attenuation values at this energy
	mu_comp_at_E = zeros(1, noMaterials);
	mu_by_rho_comp_at_E = zeros(1, noMaterials);
	mu_by_rho_en_comp_at_E = zeros(1, noMaterials);
	for m = 1:noMaterials
		mu_comp_at_E(m) = interp1geom(E_comp{m}, mu_comp{m}, E_alloy(e));
		mu_by_rho_comp_at_E(m) = interp1geom(E_comp{m}, mu_by_rho_comp{m}, E_alloy(e));
		mu_by_rho_en_comp_at_E(m) = interp1geom(E_comp{m}, mu_by_rho_en_comp{m}, E_alloy(e));
	end
	
	% compute alloy's attenuation coefficient as volume-weighted attenuation coefficient (since volume translates to thickness if the area is fixed)
	mu_alloy(e) = dot(mu_comp_at_E, V_comp/sum(V_comp));
	
	% compute alloy's mass attenuation coefficient as mass-weighted mass attenuation coefficient (see NIST documentation or derive from above formula)
	mu_by_rho_alloy(e) = dot(mu_by_rho_comp_at_E, m_comp/sum(m_comp));
	
	% the mass-weighted mass attenuation coefficient should be the same as the attenuation coefficient divided by the density
	if abs(mu_by_rho_alloy(e) / (mu_alloy(e)/rho_alloy) - 1) > 10*eps
		error('%f ~= %f !', mu_by_rho_alloy(e), my_alloy(e)/rho_alloy);
	end
	
	% compute alloy's mass energy-absorption coefficient as mass-weighted mass energy-absorption coefficient
	mu_by_rho_en_alloy(e) = dot(mu_by_rho_en_comp_at_E, m_comp/sum(m_comp));
end


%% Plots (for Verification)

figure;
hold all;
for m = 1:noMaterials
	plot(E_comp{m}, mu_comp{m});
end
plot(E_alloy, mu_alloy, 'k');
xlabel('E (keV)');
ylabel('\mu (1/cm)');
title('Attenuation Coefficient');
legend({materialNames{:}, alloyName});

figure;
hold all;
for m = 1:noMaterials
	plot(E_comp{m}, mu_by_rho_comp{m});
end
plot(E_alloy, mu_by_rho_alloy, 'k');
xlabel('E (keV)');
ylabel('\mu/\rho (cm^2/g)');
title('Mass Attenuation Coefficient');
legend({materialNames{:}, alloyName});


%% Store Computed Values

filename = [materialsDir 'density_' alloyName '.txt'];
fileId = fopen(filename, 'w');
fprintf(fileId, '%.3e', rho_alloy);
fclose(fileId);

filename = [materialsDir 'attenuation_' alloyName '.txt'];
fileId = fopen(filename, 'w');
for e = 1:length(E_alloy)
	fprintf(fileId, '%.5e %.3e %.3e', E_alloy(e), mu_by_rho_alloy(e), mu_by_rho_en_alloy(e));
	if e ~= length(E_alloy), fprintf(fileId, '\n'); end
end
fclose(fileId);
