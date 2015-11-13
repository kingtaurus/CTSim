function process_seq_file

dir = 'E:\Data\NasaFlame\Nov_10_2015_Study\3ppi_interface_60kV_50mA\';

seq_filename = 'air_hot_02'
    
    DetectorSize = [1024 768];
    size_per_frame = DetectorSize(1)*DetectorSize(2);%uint16
    
    filename = 'proj.raw';
    filterindex = 2;
    pathname = [dir seq_filename '\'];
    mkdir( pathname );
    
    seq_file = [dir seq_filename '.seq'];
    
    raw_file_base = [pathname 'proj_'];
    index = 5;
    
    ReadSeqFile_raw(seq_file, raw_file_base);
    
    fprintf( '%s Done! \n' );
    
end

