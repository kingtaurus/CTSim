function phantomToGeant

% get data
[data meta] = readMetaImage;
if isequal(data, 0) && isequal(meta, 0), error('No input file specified!'); end

% transpose data to write the DICOM file with the correct orientation
data = data';
% meta.ElementSpacing = meta.ElementSpacing([2 1]);

% get other parameters
sizeImg = size(data);
if length(sizeImg) == 2
	sizeImg(3) = 1;
end
spacingImg = meta.ElementSpacing;
if length(spacingImg) == 2
	spacingZ = input('What is the spacing in z direction (in mm)? ');
	if isempty(spacingZ) || ~isscalar(spacingZ), error('Input is not a scalar!'); end
	spacingImg(3) = spacingZ;
end
offsetZ = input('What is the offset of the volume''s center from the coordinate\nsystem''s center in z direction (in mm, default = 0mm)? ');
if isempty(offsetZ), offsetZ = 0; end
if ~isscalar(offsetZ) || ~isreal(offsetZ), error('Input is not a scalar!'); end
origin = -(sizeImg-1).*spacingImg/2;
origin(3) = origin(3) + offsetZ;

% scale intensities to prepare for DICOM
data = double(data);
intensityMaxIn = max(data(:));
intensityMinIn = min(data(:));
intensityMinOut = 0;
intensityIncrementOut = 100;
intensityMaxOut = (intensityMaxIn-intensityMinIn)*intensityIncrementOut + intensityMinOut;
data = (data-intensityMinIn)/(intensityMaxIn-intensityMinIn); % map to [0, 1]
data = data*(intensityMaxOut-intensityMinOut) + intensityMinOut; % map to [minOut maxOut]
data = data/(2^16-1); % convert to DICOM range

% ask for DICOM filename
filters = { ...
	'*.dcm',  'DICOM single file (*.dcm)'; ...
	'*.*',  'All Files (*.*)' ...
	};
[file, path] = uiputfile(filters, 'Save DICOM Image');
if isscalar(file) && isscalar(path) && isequal(file, 0) && isequal(path, 0)
	disp('No filename specified - aborting.');
	return;
else
	filenameDicom = [path, file];
end

% preliminary write to create meta information structure
dicomwrite(data, filenameDicom);
info = dicominfo(filenameDicom);

% augment DICOM meta information
info.Modality = 'CT';
info.SliceThickness = spacingImg(3);
info.SpacingBetweenSlices = spacingImg(3);
% info.ReconstructionDiameter = 
info = rmfield(info, 'PatientOrientation');
info.ImagePositionPatient = origin;
info.ImageOrientationPatient = [1 0 0 0 1 0];
info.SliceLocation = -info.ImagePositionPatient(3);
info.PixelSpacing = spacingImg(1:2)';
info.RescaleIntercept = 0;
info.RescaleSlope = 1;

% final write with complete DICOM meta information
dicomwrite(data, filenameDicom, info, 'CreateMode', 'copy');

% print list of DICOM values
data = dicomread(filenameDicom);
disp('List of DICOM intensities:');
disp(unique(data(:)));
