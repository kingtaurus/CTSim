close all;
clear all;

[filename, dir, filterindex] = uigetfile( ...
    {  '*.seq','Sequence-files (*.seq)'; ...
    '*.viv','Individual viva file (*.viv)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Pick a sequence file you want to separate: ', ...
    'MultiSelect', 'on');

if dir == 0
    return;
end;

DetectorSize = [1024 768];
size_per_frame = DetectorSize(1)*DetectorSize(2);%uint16
seq_file = strcat(dir,'', filename);

[filename, pathname, filterindex] = uiputfile( ...
    {'*.viv', 'Individual viva file (*.viv)';...
    '*.raw','Raw data file (*.raw)';...
    '*.*',  'All Files (*.*)'},...
    'Save as separate viv file (base name required):');

if isequal(filename,0) || isequal(pathname,0)
    disp('User selected Cancel');
    return;
end;

if filterindex == 1
    index = strfind(filename, '.viv');
    viv_file_base = strcat(pathname, filename(1:index-1), '_');
    ReadSeqFile_viv(seq_file, viv_file_base);
else if filterindex == 2
        index = strfind(filename, '.raw');
        raw_file_base = strcat(pathname, filename(1:index-1), '_');
        ReadSeqFile_raw(seq_file, raw_file_base);
    end;
end;

