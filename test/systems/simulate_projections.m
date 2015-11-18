load 'temp.mat';

% Load simulation parameters and datas
[phan, map] = loadXCATPhantom(p);

[ geom ] = loadProjectionGeometryCT( p );
spectrum = loadSpectraCT(p, geom, 2e6);
sinoPhotonCounts= simulatePhotonCountingData( phan, geom, spectrum, sinosDirKeV, false );


turns = 4;
pitch = 0.45;
[ geom ] = loadProjectionGeometryHelicalCT( p, turns, pitch );
sinosHelicalDirKeV = [sinosDirKeV sprintf( '%i-turns-%i-pitch', turns, round( pitch * 100 ) ), '/'];
spectrum = loadSpectraCT(p, geom, 2e5);
sinoPhotonCounts= simulatePhotonCountingData( phan, geom, spectrum, sinosHelicalDirKeV, false  );


turns = 4;
pitch = 0.95;
[ geom ] = loadProjectionGeometryHelicalCT( p, turns, pitch );
sinosHelicalDirKeV = [sinosDirKeV sprintf( '%i-turns-%i-pitch', turns, round( pitch * 100 ) ), '/'];
spectrum = loadSpectraCT(p, geom, 2e5);
sinoPhotonCounts= simulatePhotonCountingData( phan, geom, spectrum, sinosHelicalDirKeV, false  );