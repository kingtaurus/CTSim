dir = 'E:\Data\NasaFlame\Nov_5_2015_Study\';
seq_filename = 'DiffusionFlameLaminar_1';


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

fprintf( 'Done!' );