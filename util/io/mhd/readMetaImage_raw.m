function data = readMetaImage_raw(filenameRaw, DimSize, ElementType, byteOrderMSB)
% data = readMetaImage_raw(filenameRaw, DimSize, ElementType, byteOrderMSB)
%   is a helper function and reads the RAW data file using the given
%   parameters.
%
% Instead of using this function, better use readMetaImage for
% reading your files.
%
% Andreas Keil, keila@in.tum.de, 2008-12-17
% Modified for .mha format : Mehmet Yigitsoy, yigitsoy@in.tum.de, 2010-06-01


%% Check input

% parse file name
if ~exist(filenameRaw, 'file')
	error('The given raw file does not exist!');
end

% parse DimSize
if ~(isvector(DimSize) && min(DimSize) > 0)
	error('The given DimSize is not valid!');
end

nofElem = prod(DimSize);

% parse ElementType
switch ElementType
	case 'MET_UCHAR'
		datatype_matlab = 'uint8';
        offset = -nofElem;
	case 'MET_CHAR';
		datatype_matlab = 'int8';
        offset = -nofElem;
	case 'MET_SHORT';
		datatype_matlab = 'int16';
        offset = -nofElem*2;
	case 'MET_USHORT';
		datatype_matlab = 'uint16';
        offset = -nofElem*2;
	case 'MET_FLOAT';
		datatype_matlab = 'single';
        offset = -nofElem*4;
	case 'MET_DOUBLE';
		datatype_matlab = 'double';
        offset = -nofElem*8;
	otherwise
		error('Unknown ElementType!');
end

% parse byteOrderMSB
if nargin < 4
	byteOrderMSB = false;
elseif ~(isscalar(byteOrderMSB) && islogical(byteOrderMSB))
	error('Specified byteOrderMSB is not valid!');
end
if byteOrderMSB
	error('Specified byteOrderMSB is not implemented!');
end

%% Read raw file
fileId = fopen(filenameRaw);
fseek(fileId,offset,'eof'); % in order the skip the text part of mha file
data = fread(fileId, nofElem, [datatype_matlab '=>' datatype_matlab]);
if length(DimSize) > 1
	data = reshape(data, DimSize);
end
fclose(fileId);
