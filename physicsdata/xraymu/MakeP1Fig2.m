% makeP1Fig2.m
% Make Fig. 2 of Post 1
% REA 2011-05-28

%% make  Figure
xray_energies = 10:0.2:10000; % keV 
mus = []; % will hold mu data for xray_energies along rows and for all elements in material along columns 
  %plot the results
fh = figure;  
dat = BodyMaterialCompositionFunc; % composition of body materials
[weights,chemical_symbols] = XrayData;
for name = {'bone_cortical'}  % can add other tissues
  dz = dat.(char(name)); % table with composition and fraction by weight
  mu_material = xraymu(dz(2:end,:),xray_energies);
  for kelem = 2:size(dz,1)
    mu_z = xraymu([dz(kelem,1) 1],xray_energies);
    mus = [mus mu_z];
  end
  loglog(xray_energies,mus,'-k');
  hold on
  [tmp,idx] = min(abs(xray_energies-15));
  xtext = 1.05*xray_energies(idx);
  for kelem = 2:size(dz,1)
    ytext = mus(idx,kelem-1);
    symbol = chemical_symbols(dz(kelem,1),:);
       % S and P fall on top of each other,so ...
    if strcmp(symbol(1),'S')
       ytext = 1.1*ytext;
    end
    if strcmp(symbol(1),'P')
       ytext = 0.9*ytext;
    end
    ht = text(xtext,ytext,symbol,'fontsize',12);
  end
  hold off
end
set(gca,'xlim',[10 10000],'fontsize',12);
set(fh,'Color','white','Position',[151   937   565   532],'name',mfilename);
xlabel('X-ray Energy (keV)');
ylabel('Mass Attenuation Coefficient (cm^2/gm)');


return;
%% export it as eps for pdf and png for html
exportfig(fh,'E:\Projects\Blog\Posts\P2xraymu\MusElemsInBody.eps','width',4,'LockAxes',0);
figno = sprintf('-f%d',fh);
print(figno,'E:\Projects\Blog\Posts\P2xraymu\MusElemsInBody.png','-dpng','-cmyk') 