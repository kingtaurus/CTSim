function writeMetaImage(data, varargin)
% writeMetaImage  writes data in the MetaImage (MHD+RAW) format.
%
% writeMetaImage(data, filename)
%   or
% writeMetaImage(data, filename, meta)
%   or
% writeMetaImage(data, filename, meta, forceMetaProperties)
%   or
% writeMetaImage(data, filename, sp)
%   or
% writeMetaImage(data, filename, sp, orig)
%   takes a data image or volume and a filename (with or without path
%   and with ".mhd" or without extension) and writes the .mhd/.raw
%   file pair, using optionally given meta information (see
%   writeMetaImage_mhd for details of the meta structure).
%
% If no filename is given, writeMetaImage switches to interactive mode and
% displays a file save dialog for selecting the MHD file to write to.
%
% Examples:
%   writeMetaImage(data); % interactive mode for the file name, no meta data
%   writeMetaImage(data, meta); % interactive mode for the file name
%   writeMetaImage(data, 'head.mhd'); % file name given, but no spacing or origin
%   writeMetaImage(data, 'somedir/head.mhd', meta);
%   writeMetaImage(data, 'head.mhd', spacing); % no origin information written
%   writeMetaImage(data, 'head.mhd', spacing, origin);
%
% Attention: If you specify a raw file name prefix in the meta data which
%   does not correspond to the mhd file's prefix (e.g. when reading a
%   MetaImage and writing it to another file using the same meta data
%   afterwards), there is a conflict which you have to resolve using the
%   boolean parameter "forceMetaProperties". If you set this to "true", the
%   raw file name stored in the meta data will be used (which may or may
%   not correspond to the mhd file name in this case). If you specify
%   "false" (the default value), a raw file name corresponding to the mhd
%   file name will be used.
%   The same goes for the data type. Setting "forceMetaProperties" to
%   "true" results in a data conversion if necessary.
%
% See also readMetaImage.
%
% Uses writeMetaImage_mhd and writeMetaImage_raw.
%
% Andreas Keil, keila@in.tum.de, 2008-12-17


%% Input check
if nargin < 1
	error('You have to specifiy at least the data variable for writing a MetaImage!');
end
if nargin < 2 || ~ischar(varargin{1})
	% file name not specified, display dialog box
	filters = { ...
		'*.mhd',  'MetaImage (*.mhd)'; ...
		'*.*',  'All Files (*.*)' ...
		};
	[file, path] = uiputfile(filters, 'Save MetaImage');
	if isscalar(file) && isscalar(path) && isequal(file, 0) && isequal(path, 0)
		disp('No filename specified - aborting.');
		return;
	else
		% add selected filename to argument list
		varargin = {[path, file], varargin{:}};
		argsGiven = nargin + 1;
	end
else
	argsGiven = nargin;
end


%% Parse standard input

% deduce dimensionality and size from data (choosing the minimum dimensionality as temporary guess)
if isscalar(data)
	meta_tmp.NDims = 1;
	meta_tmp.DimSize = 1;
elseif isvector(data)
	meta_tmp.NDims = 1;
	meta_tmp.DimSize = length(data);
else
	meta_tmp.NDims = ndims(data);
	meta_tmp.DimSize = size(data);
end

% parse filename
[path, name, ext] = fileparts(varargin{1});
if ~isempty(path)
	path = [path '/'];
end
if ~(strcmpi(ext, '.mhd') || isempty(ext))
	error('The given extension is not supported!');
end
filenameRaw = [name '.raw'];
filenameMhd = [name '.mhd'];
clear name ext; % keep only path, filenameMhd, and filenameRaw
meta_tmp.ElementDataFile = filenameRaw;

% parse data type
switch class(data)
	case 'uint8'
		meta_tmp.ElementType = 'MET_UCHAR';
	case 'int8'
		meta_tmp.ElementType = 'MET_CHAR';
	case 'int16'
		meta_tmp.ElementType = 'MET_SHORT';
	case 'uint16'
		meta_tmp.ElementType = 'MET_USHORT';
	case 'single'
		meta_tmp.ElementType = 'MET_FLOAT';
	case 'double'
		meta_tmp.ElementType = 'MET_DOUBLE';
	otherwise
		% for different possible data types see help fread
		warning('Datatype unknown - converting to float!');
		meta_tmp.ElementType = 'MET_FLOAT';
		data = cast(data, 'single');
end


%% Parse further input
if argsGiven >= 3 && isstruct(varargin{2})
	% parse option for raw file name specification
	if argsGiven >= 4
		if ~islogical(varargin{3}) || ~isscalar(varargin{3})
			error('The "forceMetaProperties" argument has to be a single boolean!');
		end
		forceMetaProperties = varargin{3};
	else
		forceMetaProperties = false;
	end
	
	% parse given meta data structure
	meta = varargin{2};

	if isfield(meta, 'NDims')
		if meta.NDims < meta_tmp.NDims
			error('Given NDims is smaller than the given data''s dimensionality!');
		elseif meta.NDims > meta_tmp.NDims
			meta_tmp.DimSize = [meta_tmp.DimSize ones(1, meta.NDims-meta_tmp.NDims)]; % append singleton dimensions to temporarily guessed size
		end
	else
		meta.NDims = meta_tmp.NDims;
	end

	if isfield(meta, 'DimSize')
		if ~isequal(meta.DimSize, meta_tmp.DimSize) 
			error('Given DimSize does not match data''s size!');
		end
	else
		meta.DimSize = meta_tmp.DimSize;
	end

	if isfield(meta, 'ElementType') && forceMetaProperties
		% use the given data type and convert data accordingly
		data = cast(data, meta.ElementType);
	else
		% use the computed data type
		meta.ElementType = meta_tmp.ElementType;
	end

	if isfield(meta, 'BinaryDataByteOrderMSB') && meta.BinaryDataByteOrderMSB
		error('Given BinaryDataByteOrderMSB is not implemented!');
	end

	if isfield(meta, 'ElementByteOrderMSB') && meta.ElementByteOrderMSB
		error('Given ElementByteOrderMSB is not implemented!');
	end

	if isfield(meta, 'ElementDataFile') && forceMetaProperties
		% use the given file name
	else
		% use the computed file name
		meta.ElementDataFile = meta_tmp.ElementDataFile; % use computed raw file name
	end
	
else
	% parse optionally given spacing and origin vectors
	meta = meta_tmp;

	if argsGiven >= 3
		if ~isvector(varargin{2}) || length(varargin{2}) < meta.NDims
			error('The "sp" argument has to be at least a %d-vector of spacings!', meta.NDims);
		end
		meta.ElementSpacing = varargin{2};
		if length(varargin{2}) > meta.NDims
			meta.NDims = length(varargin{2});
			meta.DimSize = [meta.DimSize ones(1, meta.NDims-length(meta.DimSize))]; % append singleton dimensions to size
		end
	end

	if argsGiven >= 4
		if ~isvector(varargin{3}) || (length(varargin{3}) ~= meta.NDims)
			error('The "orig" argument has to be a %d-vector of position values!', meta.NDims);
		end
		meta.Position = varargin{3};
	end

end


%% Write mhd file
writeMetaImage_mhd([path filenameMhd], meta);


%% Write raw file
writeMetaImage_raw([path meta.ElementDataFile], data);
