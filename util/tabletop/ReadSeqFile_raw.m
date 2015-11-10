function ReadSeqFile_raw(fname_seq, fname_raw_base, type, bTranspose)
%ReadSeqFile       Read .seq file
%   data = ReadSeqFile(fname, type, bTranspose = 1)
%
%   type : either 'int32' so output data type is int32, or default in which
%   data is of the double type
%
%   Note: data = permute(data_read_out, [1 3 2]);
%   thus data is 3D and arranged in (z,y,x)
%
%   $Id: ReadSeqFile.m,v 1.4 2009/07/30 20:48:35 msun Exp $


% .viv (and .seq) header is 2048 bytes long
VIV_HEADERSIZE = 2048 / 4;

fid = fopen(fname_seq);
if fid < 0
    error([fname_seq ' does not exist.'])
end

% read header
header = fread(fid, VIV_HEADERSIZE, 'int32');

% get the image size
nx = header(5);
ny = header(6);
nz = header(8);

data = zeros(1,nx,ny,'int32');

angle_inc = rad2deg(2*pi/nz);
% read the data
for n = 1 : nz
    data(1,:,:) = fread(fid, [nx, ny], 'ushort=>int32');

    % decode if file is in 2MSB format
    dataType = header(460);
    if bitand(dataType,512),
        index1 = find(data >= 49152);
        index2 = find(data >= 32768 & data < 49152);
        index3 = find(data >= 16384 & data < 32768);
        data(index1) = (data(index1) - 49152) * 8;
        data(index2) = (data(index2) - 32768) * 4;
        data(index3) = (data(index3) - 16384) * 2;
    end

    % convert to double by default
    if ~(exist('type','var') && strcmp(type,'int32')),
        data = double(data);
    end

    % transpose to fit into Matlab's convention
    if ~(exist('bTranspose','var') && bTranspose==0),
        data = permute(data, [1 3 2]);
    end

    
    frame = reshape(data(1,:,:), [ny nx]);
   
    raw_file_name = strcat(fname_raw_base, num2str(n-1, '%05d'), '.raw');
    
    fid_raw = fopen(raw_file_name, 'w');
    fwrite(fid_raw, frame', 'float32');
    fclose(fid_raw);
        
    clear data;
end
fclose(fid);