function imgPhantom = generatePhantomHeadMaterials2Patterns2Fillings(overSampling, outputFilename)


%% Parse Input Parameters

if nargin < 2
	% ask for filename
	filters = { ...
		'*.mhd',  'MetaImage files (*.mhd)'; ...
		'*.*',  'All Files (*.*)' ...
		};
	[file, path] = uiputfile(filters, 'Save Phantom Image');
	if isscalar(file) && isscalar(path) && isequal(file, 0) && isequal(path, 0)
		disp('No filename specified - aborting.');
		return;
	else
		outputFilename = [path, file];
	end
end

if nargin < 1
	overSampling =1;
end


%% Define Ellipses

% columns 1,2: ellipse's center (mm)
% columns 3,4: ellipse's axes (mm)
% column  5:   ellipse's rotation (degrees)
% column  6:   material index (in powers of 2, since the "ellipses" method accumulates values)

ellipsesParams = [
	% head
	  0   0   80   88      0  1;
	% fillings
	-44 -28.8  1.6  1.6    0  2;
	 36 -48    1.6  1.6    0 16;
	% line pattern between fillings
	  0 -44.8  9.6  0.8  -15  4;
	  0 -42.4  9.6  0.8  -15  4;
	  0 -40.0  9.6  0.8  -15  4;
	  0 -37.6  9.6  0.8  -15  4;
	  0 -35.2  9.6  0.8  -15  4;
	% teeth
	-44 -28.8  9.6  4.8  -75  8;
	-36 -48    9.6  4.8  -50  8; 
	-20 -60    9.6  4.8  -30  8; 
	  0 -64    9.6  4.8    0  8; 
	 20 -60    9.6  4.8   30  8; 
	 36 -48    9.6  4.8   50  8; 
	 44 -28.8  9.6  4.8   75  8; 
	% jaw
	-52  20   16    3.2 -100  8;
	 52  20   16    3.2  100  8;
	% spine
	  0  48    9.6 9.6     0  8;
	% line pattern away from fillings
	  0   3.2  9.6 0.8   -15  4;
	  0   5.6  9.6 0.8   -15  4;
	  0   8    9.6 0.8   -15  4;
	  0  10.4  9.6 0.8   -15  4;
	  0  12.8  9.6 0.8   -15  4;
];
% negate y coordinates and rotation angles (due to weird convention, the
% "ellipses" method draws on an inverted y axis, but we wanted to specify
% the coordinates and angles in the mathematical convention above)
ellipsesParams(:, [2 5]) = -ellipsesParams(:, [2 5]);


%% Define Phantom Size and Spacing

nx = round(280 * overSampling); % size in x direction (pixels)
ny = round(300 * overSampling); % size in y direction (pixels)
dx = 0.6 / overSampling; % spacing in x direction (mm)
dy = 0.6 / overSampling; % spacing in y direction (mm)


%% Generate Phantom from Ellipses Definitions

imgPhantom = ellipses(nx, ny, ellipsesParams, dx, dy);


%% Re-label Material Indexes to Consecutive Numbers

imgPhantom(imgPhantom == 1) = 2;
imgPhantom(imgPhantom == 5) = 1;
imgPhantom(imgPhantom == 9) = 3;
imgPhantom(imgPhantom == 11) = 5;
imgPhantom(imgPhantom == 25) = 4;


%% Save Phantom Image File

writeMetaImage(imgPhantom, outputFilename, [dx dy]);
