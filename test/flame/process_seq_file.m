function process_seq_file( dir, seq_filename )



DetectorSize = [1024 768];
size_per_frame = DetectorSize(1)*DetectorSize(2);%uint16

filename = 'proj.raw';
filterindex = 2;

pathname = [dir seq_filename '\'];

if exist( pathname, 'dir' ) 
    fprintf('Dirtory %s  already exist. \n', pathname);
    return;
end

mkdir( pathname );

seq_file = [dir seq_filename '.seq'];

raw_file_base = [pathname 'proj_'];

ReadSeqFile_raw(seq_file, raw_file_base);

fprintf( '%s Done! \n' );

end

