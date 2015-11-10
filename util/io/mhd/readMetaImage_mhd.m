function meta = readMetaImage_mhd(filenameMhd);
% readMetaImage_mhd(filenameMhd)  is a helper function and
%   reads the MHD header file using the given file name and meta
%   information.
%
% Instead of using this function, better use readMetaImage for
% reading your files.
%
% Andreas Keil, keila@in.tum.de, 2008-12-17


%% Input check
[path, name, ext] = fileparts(filenameMhd);
clear path name; % keep only ext
if ~strcmpi(ext, '.mhd')
	error('The given extension (of file "%s") is not supported!', filenameMhd);
end
if ~exist(filenameMhd, 'file')
	error('The given mhd file "%s" does not exist!', filenameMhd);
end


%% Open file
fileId = fopen(filenameMhd, 'r');


%% Parse file
while true

	% read a text line
	textline = fgetl(fileId);
	if ~ischar(textline)
		break;
	end

	% make sure the meta data structure exists, so that fields can also be added using setfield()
	meta.dummy = 1;
	meta = rmfield(meta, 'dummy');
	
	% parse the line
	[name, textline] = strtok(textline);
	[equalsign, values] = strtok(textline);
	if isempty(name) || isempty(equalsign) || ~strcmpi(equalsign, '=') || isempty(values);
		error('File parsing error!');
	end

	switch name
		case 'NDims'
			parsed = textscan(values, '%d');
			if length(parsed) ~= 1, error('Parsing error!'); end
			meta = setfield(meta, 'NDims', double(parsed{1}));
			if ~(isscalar(meta.NDims) && length(meta.NDims) == 1 && meta.NDims > 0)
				error('Illegal NDims!');
			end
		case 'DimSize'
			parsed = textscan(values, '%d', meta.NDims);
			if length(parsed) ~= 1, error('Parsing error!'); end
			meta.DimSize = double(parsed{1})';
			if ~(isvector(meta.DimSize) && length(meta.DimSize) == meta.NDims && min(meta.DimSize) > 0)
				error('Illegal DimSize!');
			end
		case 'BinaryDataByteOrderMSB'
			parsed = textscan(values, '%s');
			if length(parsed) ~= 1, error('Parsing error!'); end
			meta.BinaryDataByteOrderMSB = parsed{1}{1};
			switch meta.BinaryDataByteOrderMSB
				case 'False'
					meta.BinaryDataByteOrderMSB = false;
				case 'True'
					meta.BinaryDataByteOrderMSB = true;
				otherwise
					error('Illegal BinaryDataByteOrderMSB!');
			end
		case 'ElementSpacing'
			parsed = textscan(values, '%f', meta.NDims);
			if length(parsed) ~= 1, error('Parsing error!'); end
			meta.ElementSpacing = double(parsed{1})';
			if ~(isvector(meta.ElementSpacing) && length(meta.ElementSpacing) == meta.NDims && min(meta.ElementSpacing) >= 0)
				error('Illegal ElementSpacing!');
			end
			if sum(meta.ElementSpacing == 0) > 0
				warning('Zero ElementSpacing component.');
			end
		case 'Position'
			parsed = textscan(values, '%f', meta.NDims);
			if length(parsed) ~= 1, error('Parsing error!'); end
			meta.Position = double(parsed{1})';
			if ~(isvector(meta.Position) && length(meta.Position) == meta.NDims)
				error('Illegal Position!');
			end
		case 'ElementType'
			parsed = textscan(values, '%s');
			if length(parsed) ~= 1, error('Parsing error!'); end
			meta.ElementType = parsed{1}{1};
			if ~(ischar(meta.ElementType) && strncmp(meta.ElementType, 'MET_', 4))
				error('Illegal ElementType!');
			end
		case 'ElementByteOrderMSB'
			parsed = textscan(values, '%s');
			if length(parsed) ~= 1, error('Parsing error!'); end
			meta.ElementByteOrderMSB = parsed{1}{1};
			switch meta.ElementByteOrderMSB
				case 'False'
					meta.ElementByteOrderMSB = false;
				case 'True'
					meta.ElementByteOrderMSB = true;
				otherwise
					error('Illegal ElementByteOrderMSB!');
			end
		case 'ElementDataFile'
			parsed = textscan(values, '%s');
			if length(parsed) ~= 1, error('Parsing error!'); end
			meta.ElementDataFile = parsed{1}{1};
			if ~ischar(meta.ElementDataFile)
				error('Illegal ElementDataFile!');
			end
		otherwise
			values = strtrim(values);
 			meta = setfield(meta, name, values);
 			if ~ischar(getfield(meta, name))
 				error('Parsing Error!');
 			end
			%disp(sprintf('Storing unknown tag "%s = %s".', name, getfield(meta, name)));
	end % switch
	
end % while true


%% Close file
fclose(fileId);
