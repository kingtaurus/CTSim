function writeMetaImage_mhd(filenameMhd, meta)
% writeMetaImage_mhd(filenameMhd, meta)  is a helper function and
%   writes the MHD header file using the given file name and meta
%   information.
%
% Instead of using this function, better use writeMetaImage for
% writing your files.
%
% Andreas Keil, keila@in.tum.de, 2008-12-17


%% parse file name
[path, name, ext] = fileparts(filenameMhd);
if ~isempty(path)
	path = [path '/'];
end
if ~(strcmpi(ext, '.mhd') || isempty(ext))
	error('The given extension is not supported!');
end
filenameMhd = [path name '.mhd'];

%% Open file
mhdFile = fopen(filenameMhd, 'w');


%% Write meta data (Attention: The order in which the meta data is written matters!)

% store number of dimensions for further tests
if ~(isfield(meta, 'NDims') && isscalar(meta.NDims) && meta.NDims >= 1)
	error('You have to specify a correct NDims value!');
end
N = meta.NDims;

% NDims is mandatory
if ~(isfield(meta, 'NDims') && isscalar(meta.NDims) && meta.NDims >= 1)
	error('You have to specify a correct NDims value!');
end
fprintf(mhdFile, 'NDims = %d\n', meta.NDims);
meta = rmfield(meta, 'NDims');

% DimSize is mandatory
if ~(isfield(meta, 'DimSize') && isvector(meta.DimSize) && length(meta.DimSize) == N && min(meta.DimSize) > 0)
	error('You have to specify a correct DimSize of the data!');
end
fprintf(mhdFile, 'DimSize =');
for dim = 1:N
	fprintf(mhdFile, ' %d', meta.DimSize(dim));
end
fprintf(mhdFile, '\n');
meta = rmfield(meta, 'DimSize');

% byte order is optional
if isfield(meta, 'BinaryDataByteOrderMSB')
	if ~(isscalar(meta.BinaryDataByteOrderMSB) && islogical(meta.BinaryDataByteOrderMSB))
		error('The given BinaryDataByteOrderMSB is not valid!');
	end
	if meta.BinaryDataByteOrderMSB
		fprintf(mhdFile, 'BinaryDataByteOrderMSB = True\n');
	else
		fprintf(mhdFile, 'BinaryDataByteOrderMSB = False\n');
	end
	meta = rmfield(meta, 'BinaryDataByteOrderMSB');
end

% spacing is mandatory and defaults to 1
if ~isfield(meta, 'ElementSpacing')
	meta.ElementSpacing = ones(N, 1);
end
if ~(isvector(meta.ElementSpacing) && length(meta.ElementSpacing) == N && min(meta.ElementSpacing) >= 0)
	error('You have to specify a correct ElementSpacing of the data!');
end
fprintf(mhdFile, 'ElementSpacing =');
for dim = 1:N
	fprintf(mhdFile, ' %f', meta.ElementSpacing(dim));
end
fprintf(mhdFile, '\n');
meta = rmfield(meta, 'ElementSpacing');

% position is optional
if isfield(meta, 'Position')
	if ~(isvector(meta.Position) && length(meta.Position) == N)
		error('You have to specify a correct Position of the data!');
	end
	fprintf(mhdFile, 'Position =');
	for dim = 1:N
		fprintf(mhdFile, ' %f', meta.Position(dim));
	end
	fprintf(mhdFile, '\n');
	meta = rmfield(meta, 'Position');
end

% data type is mandatory
if ~(isfield(meta, 'ElementType') && ischar(meta.ElementType) && strcmp(meta.ElementType(1:4), 'MET_'))
	error('You have to specify a correct ElementType!');
end
fprintf(mhdFile, 'ElementType = %s\n', meta.ElementType);
meta = rmfield(meta, 'ElementType');

% byte order is optional
if isfield(meta, 'ElementByteOrderMSB')
	if ~(isscalar(meta.ElementByteOrderMSB) && islogical(meta.ElementByteOrderMSB))
		error('The given ElementByteOrderMSB is not valid!');
	end
	if meta.ElementByteOrderMSB
		fprintf(mhdFile, 'ElementByteOrderMSB = True\n');
	else
		fprintf(mhdFile, 'ElementByteOrderMSB = False\n');
	end
	meta = rmfield(meta, 'ElementByteOrderMSB');
end

% raw file name is mandatory
if ~(isfield(meta, 'ElementDataFile') && ischar(meta.ElementDataFile) && strcmpi(meta.ElementDataFile(end-3:end), '.raw'))
	error('You have to specify a correct ElementDataFile!');
end
fprintf(mhdFile, 'ElementDataFile = %s\n', meta.ElementDataFile);
meta = rmfield(meta, 'ElementDataFile');

% write other unknown fields
names = fieldnames(meta);
for i = 1:length(names)
	field = names(i);
	field = field{1};
	value = getfield(meta, field);
	disp(sprintf('Writing unknown tag "%s = %s".', field, value));
	fprintf(mhdFile, '%s = %s\n', field, value);
end

%% Close file
fclose(mhdFile);
