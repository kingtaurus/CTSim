function writeMetaImage_raw(filenameRaw, data)
% writeMetaImage_raw(filenameRaw, data)  is a helper function and
%   writes the RAW data file using the given file name and data.
%
% Instead of using this function, better use writeMetaImage for
% writing your files.
%
% Andreas Keil, keila@in.tum.de, 2008-12-17


%% Parse file name
[path, name, ext] = fileparts(filenameRaw);
if ~isempty(path)
	path = [path '/'];
end
if ~(strcmpi(ext, '.raw') || isempty(ext))
	error('The given extension is not supported!');
end
filenameRaw = [path name '.raw'];

%% Open file
rawFile = fopen(filenameRaw, 'w');

%% Write data
fwrite(rawFile, data(:), class(data));

%% Close file
fclose(rawFile);
