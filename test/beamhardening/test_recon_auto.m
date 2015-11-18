% load reconstruction parameters
load 'temp.mat';

geom = loadProjectionGeometryCT( p );

spectrum = loadSpectraCT(p, geom, 2e6);

dir = 'E:\Data\NasaFlame\Nov_10_2015_Study\10ppi_interface_60kV_50mA\';

file_list = { 'air_01'; 'air_02'; 'background_01'; ...
    'background_02'; 'burn_01'; 'burn_02'; 'burn_03'; 'burn_04'; 'burn_05'; 'kr100_01'; 'kr100_02' };

for i = 1 : length( file_list )
    
    seq_filename = file_list{i};
    
    dataPath = [dir seq_filename '\'];
    
    process_seq_file( dir, seq_filename )
    
    sinoAtt = loadTableTopData( dataPath, geom );
    
    sinoAtt = beamHardeningMaterialCorrection(sinoAtt, spectrum, 'Quartz', 10 );
    
    img = reconFBP( sinoAtt, geom, 'ram-lak' );
    
    save([seq_filename '.mat'], 'img' );
    
end



