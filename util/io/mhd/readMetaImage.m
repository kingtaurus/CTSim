function [data, meta] = readMetaImage(filename)
% readMetaImage  reads data in the MetaImage (MHD+RAW or MHA) format.
%
% data = readMetaImage(filename)
%   or
% [data, meta] = readMetaImage(filename)
%   takes a filename and reads .mha file or the .mhd/.raw file pair. The filename
%   can be given with ".mhd",".mha" or without extension. The output
%   variable meta is optionally filled with the meta data available in
%   the meta file.
%
% If no filename is given, readMetaImage switches to interactive mode and
% displays a file open dialog for selecting the MHD or MHA files to read.
%
% Examples:
%   data = readMetaImage; % interactive mode
%   [data, meta] = readMetaImage('head.mhd');
%   data = readMetaImage('somedir/head.mhd'); % Beware: Important meta information (such as the spacing) are not available if the optional output is not used!
%
% See also writeMetaImage.
%
% Uses readMetaImage_mhd, readMetaImage_mha and readMetaImage_raw.
%
% Andreas Keil, keila@in.tum.de, 2008-12-17
% Modified : Mehmet Yigitsoy, yigitsoy@in.tum.de 2010-06-01


%% Output check
if nargout < 1
	error('You have to specify at least one output parameter for storing the data!\nE.g. use "data = readMetaImage;".');
end

%% Input check
if nargin < 1
	filters = { ...
		'*.mhd',  'MetaImage (*.mhd)'; ...
        '*.mha',  'MetaImage (*.mha)'; ...
		'*.*',  'All Files (*.*)' ...
		};
	[file, path] = uigetfile(filters, 'Open MetaImage');
	if isscalar(file) && isscalar(path) && isequal(file, 0) && isequal(path, 0)
		disp('No filename specified - aborting.');
		data = 0; meta = 0;
		return;
	else
		filename = [path, file];
	end
end


%% Parse input
[path, name, ext] = fileparts(filename);
if ~isempty(path)
	path = [path '/'];
end
if ~(strcmpi(ext, '.mhd') || strcmpi(ext, '.mha') || isempty(ext))
	error('The given extension is not supported!');
end
filenameMetaImage = [name ext];
% clear name ext; % leave around only path, filenameMhd, and filenameRaw


%% Read Meta Image File
if strcmpi(ext, '.mhd')
    % read mhd file
    meta = readMetaImage_mhd([path filenameMetaImage]);
elseif strcmpi(ext, '.mha')
	% read mha file
	meta = readMetaImage_mha([path filenameMetaImage]);
else
	error('The extension of the given file name "%s" is not supported!', filenameMetaImage);
end

%% Check meta data

% number of dimensions is mandatory
if ~(isfield(meta, 'NDims') && isscalar(meta.NDims) && meta.NDims >= 1)
	error('Specified NDims is not valid!');
end

% size is mandatory
if ~(isfield(meta, 'DimSize') && isvector(meta.DimSize) && length(meta.DimSize) == meta.NDims && min(meta.DimSize) > 0)
	error('Specified DimSize is not valid!');
end

% spacing is optional
if isfield(meta, 'ElementSpacing') && ~(isvector(meta.ElementSpacing) && length(meta.ElementSpacing) == meta.NDims && min(meta.ElementSpacing) >= 0)
	error('Specified ElementSpacing is not valid!');
end

% position is optional
if isfield(meta, 'Position') && ~(isvector(meta.Position) && length(meta.Position) == meta.NDims)
		error('Specified Position is not valid!');
end

% data type is mandatory
if ~(isfield(meta, 'ElementType') && ischar(meta.ElementType) && strcmp(meta.ElementType(1:4), 'MET_'))
	error('Specified ElementType is not valid!');
end

% byte order is optional and defaults to little endian / least significant byte first
if isfield(meta, 'ElementByteOrderMSB')
	if ~(isscalar(meta.ElementByteOrderMSB) && islogical(meta.ElementByteOrderMSB))
		error('The given ElementByteOrderMSB is not valid!');
	end
	byteOrderMSB = meta.ElementByteOrderMSB;
end
if isfield(meta, 'BinaryDataByteOrderMSB')
	if ~(isscalar(meta.BinaryDataByteOrderMSB) && islogical(meta.BinaryDataByteOrderMSB))
		error('The given BinaryDataByteOrderMSB is not valid!');
	end
	byteOrderMSB = meta.BinaryDataByteOrderMSB;
end
if isfield(meta, 'ElementByteOrderMSB') && isfield(meta, 'BinaryDataByteOrderMSB') && meta.ElementByteOrderMSB ~= meta.BinaryDataByteOrderMSB
	error('The given ElementByteOrderMSB and BinaryDataByteOrderMSB conflict!');
end
if ~isfield(meta, 'ElementByteOrderMSB') && ~isfield(meta, 'BinaryDataByteOrderMSB')
	byteOrderMSB = false;
end

% raw file name is mandatory
if ~(isfield(meta, 'ElementDataFile') && ischar(meta.ElementDataFile))
	error('The given ElementDataFile is not valid!');
else if(strcmpi(ext, '.mhd') && ~exist([path meta.ElementDataFile], 'file'))
        error('The given ElementDataFile is not valid!');
    end
end


%% Read raw file
if(strcmpi(ext, '.mhd'))
    % Read mhd file
    data = readMetaImage_raw([path meta.ElementDataFile],...
        meta.DimSize, meta.ElementType, byteOrderMSB);
else if(strcmpi(ext, '.mha'))
        % Read mha file
        data = readMetaImage_raw([path filenameMetaImage],...
            meta.DimSize, meta.ElementType, byteOrderMSB);
    end
end
