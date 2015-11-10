clear;
close all;

dataPath = '/data/headphantom-ct-nofilling-109kvp-sharp/';
info = dicominfo([ dataPath 'headphantom-ct-nofilling-000.IMA'] );

dx = info.PixelSpacing(1);
dz = info.PixelSpacing(2);
dy = info.SliceThickness;

%choose desired slices
zstart = 150;
zstop = 240;

nx = double( info.Width );
ny = double( info.Height);
nz = zstop-zstart+1;

fprintf('Read data to Matlab: \n');
BW = zeros(nx, ny, nz, 'single');

% get rough teeth
for i = zstart:zstop
    if mod(i,20) == 1
        fprintf('\tslice %d \n', i);
    end
    slice = imread([dataPath, 'headphantom-ct-nofilling.tif'], 'Index', i);
    
    % remove other bones
    I = slice > 1700;
    I(256:end,:) = 0;
    I(160:230,230:300) = 0;
    I(:,1:180) = 0;
    I(:,350:end) = 0;
    
    J = imerode( I , ones(4));
    J = imdilate( J, ones(3));
    
    slice( ~J) = 0;
    BW( :,:,i-zstart+1) = slice > 2100;
end

% display
v= BW;
[x,y,z] = meshgrid((1:nx)*dx, (1:ny)*dy, (zstart:zstop)*dz) ;
p = patch(isosurface(x,y,z,v));
isonormals(x,y,z,v,p);
set(p,'FaceColor','red','EdgeColor','none');
daspect([1,1,1])
view(3); axis tight
camlight
lighting gouraud


% write to slt
fv = isosurface(v, 0.99);
stlwrite('teeth.stl', fv);

return;
