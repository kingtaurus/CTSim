function [noTrucationKeV, noTrucationMeV] = checkFotTruncation(geomKeV, geomMeV, phan)
% Load system geometry
%
% Copyright (c) 2010-2012 by Andreas Keil, Stanford University.
% Modified by Meng Wu at 2012.9



fprintf('Check for Truncation ... \n');

noTrucationKeV = false;
noTrucationMeV = false;

widthDetectorKeV = geomKeV.detSize * geomKeV.detSpacing;
widthDetectorMeV = geomMeV.detSize * geomMeV.detSpacing;
widthProjectedLongerSideKeV = geomKeV.SDD/geomKeV.SAD * max(phan.size .* phan.spacing);
widthProjectedLongerSideMeV = geomMeV.SDD/geomMeV.SAD * max(phan.size .* phan.spacing);
widthProjectedDiagonalKeV = geomKeV.SDD/geomKeV.SAD * norm(phan.size .* phan.spacing);
widthProjectedDiagonalMeV = geomMeV.SDD/geomMeV.SAD * norm(phan.size .* phan.spacing);


if widthDetectorKeV < widthProjectedLongerSideKeV
	fprintf('keV detector does not cover projected longer side of input image. => Projections will be truncated.\n'); 
elseif widthDetectorKeV < widthProjectedDiagonalKeV
	fprintf('keV detector does not cover projected diagonal of input image. => Projections will be truncated if objects are not confined within an inscribed ellipse.\n'); 
else
    noTrucationKeV = true;
end
    
    
if widthDetectorMeV < widthProjectedLongerSideMeV
	fprintf('MeV etector does not cover projected longer side of input image. => Projections will be truncated.\n');
elseif widthDetectorMeV < widthProjectedDiagonalMeV
	fprintf('MeV detector does not cover projected diagonal of input image. => Projections will be truncated if objects are not confined within an inscribed ellipse.\n'); 
else
    noTrucationMeV = true;
end


if ( noTrucationKeV && noTrucationMeV )
    fprintf('no truncation. \n');
else
    fprintf('\n');
end
    

end



