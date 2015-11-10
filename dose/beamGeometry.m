function beamGeometry(SAD, ADD, spacingDet, noViews, fillingPositions, fillingDiameters, marginPixelsOnDetector)
% beamGeometry(SAD, ADD, spacingDet, noViews, fillingPositions, fillingDiameters, marginPixelsOnDetector)
%
% computes the geometric parameters of a set of beams that are collimated
% onto the given positions of dental fillings
%
% INPUT:
% SAD   source-to-axis distance (mm)
% ADD   axis-to-detector distance (mm)
% spacingDet   pixel spacing on the detector (mm)
% noViews   number of views / source positions to generated beams for
% fillingPositions   cell array of filling positions (cell of 2-vectors in mm)
% fillingDiameters   cell array of diameters of a fillings (mm)
% marginPixelsOnDetector   margin that is added to the left and right of a
%    fillings projection on the detector (pixels)
%
% OUTPUT:
% The user is asked for a filename to write the geometry to. If she does
% not specify a filename, the geometry parameters will be output on the
% screen instead.
% In both cases, the format will be sets of two lines per source position:
%    1 T1 N1
%    a11 c11 ... a1N c1N
%            ...
%    M TM NM
%    aM1 cM1 ... aMN cMN
% The first line specifies the view number, the angle (degrees) between the
% x-axis and the source position (denoted T... above) as well the number of
% beams for this source position (denoted N... above). The latter number is
% usually equal to the number of fillings but can be less if the beams onto
% two or more filligs overlap and therefore merge.
% The second line contains one pair of numbers per beam - the beam's center
% angle (degrees) as measured between the source-to-center line and the
% viewing direction in the mathematically positive direction (denoted a...
% above) and the beam's collimation angle (degress) which is defined as
% half the angle covered by the beam (denoted c... above).


% collimation in z
collimationHalfHeightAtIsoCenter = 5; % in mm
collimationHalfAngleZ = atan(collimationHalfHeightAtIsoCenter/SAD); % in rad

% ask for a filename to store the beam geometries in
filters = { ...
	'*.txt',  'Text files (*.txt)'; ...
	'*.*',  'All Files (*.*)' ...
	};
[file, path] = uiputfile(filters, 'Save Geometry File (cancel for display output)');
if isscalar(file) && isscalar(path) && isequal(file, 0) && isequal(path, 0)
	disp('No filename specified - writing output on screen.');
	beamOutput = 1;
else
	beamOutput = fopen([path, file], 'w');
end

% iterate over source positions
totalSteradians = 0; % cumulative beam angle (in order to be able to compute the total number of photons)
for viewNo = 1:noViews
	% source
	Theta = 2*pi/noViews * (viewNo-1);
	sourcePos = SAD * [cos(Theta); sin(Theta)];
	
	% collect opening and closing angle information for beams on all fillings
	allCollimationAngles = zeros(2*length(fillingPositions), 1);
	allCollimationAnglesType = zeros(2*length(fillingPositions), 1, 'uint8');
	for fillingNo = 1:length(fillingPositions)
		fillingPos = fillingPositions{fillingNo}(:);
		lineToFilling = fillingPos - sourcePos;
		angleToFillingWrtX = atan2(lineToFilling(2), lineToFilling(1));
		beamAngle = (pi - Theta); % default direction is to rotation center, so go first add pi to look outward, then subtract Theta to look in +x direction
		beamAngle = beamAngle + angleToFillingWrtX; % now add the angle to the filling that was compute w.r.t. the x axis
		beamAngle = beamAngle - floor(beamAngle/(2*pi))*(2*pi); % normalize to [0, 2*pi)
		fillingDiameterOnDetector = fillingDiameters{fillingNo} * (SAD+ADD) / norm(lineToFilling);
		collimationHalfWidthOnDetector = fillingDiameterOnDetector/2 + marginPixelsOnDetector*spacingDet;
		collimationHalfAngle = atan(collimationHalfWidthOnDetector/(SAD+ADD));
		allCollimationAngles(2*fillingNo-1) = beamAngle - collimationHalfAngle;
		allCollimationAnglesType(2*fillingNo-1) = 'o';
		allCollimationAngles(2*fillingNo) = beamAngle + collimationHalfAngle;
		allCollimationAnglesType(2*fillingNo) = 'c';
	end
	
	% consolidate beams for this view
	[allCollimationAngles, order] = sort(allCollimationAngles);
	allCollimationAnglesType = allCollimationAnglesType(order);
	beamAngles = [];
	collimationHalfAngles = [];
	beamsOpen = 0;
	for angleNo = 1:length(allCollimationAngles)
		if allCollimationAnglesType(angleNo) == 'o' % new beam opens
			if beamsOpen == 0 % new, separate beam -> store start angle
				angleOpen = allCollimationAngles(angleNo);
			end
			beamsOpen = beamsOpen + 1;
		elseif allCollimationAnglesType(angleNo) == 'c' % beam closes
			if beamsOpen == 1 % last open beam is closing -> store resulting beam
				angleClose = allCollimationAngles(angleNo);
				beamAngles = [beamAngles, (angleOpen+angleClose)/2]; %#ok<AGROW>
				collimationHalfAngles = [collimationHalfAngles, (angleClose-angleOpen)/2]; %#ok<AGROW>
			end
			beamsOpen = beamsOpen - 1;
		else
			error('Internal error: Unexpected collimation angle type! (Should be ''o'' for opening or ''c'' for closing anle.)');
		end
	end
	
	% write beams to geometry file
	noBeams = length(beamAngles);
	fprintf(beamOutput, '%i %f %i\n', viewNo, radtodeg(Theta), noBeams);
	for beamNo = 1:noBeams
		fprintf(beamOutput, '%f %f %f', radtodeg(beamAngles(beamNo)), radtodeg(collimationHalfAngles(beamNo)), radtodeg(collimationHalfAngleZ));
		if beamNo < noBeams
			fprintf(beamOutput, ' ');
		else
			fprintf(beamOutput, '\n');
		end
	end
	
	% accumulate total beam angle
	for beamNo = 1:noBeams
		totalSteradians = totalSteradians + 4*asin(sin(collimationHalfAngles(beamNo))*sin(collimationHalfAngleZ)); % steradian = rad^2
	end
end

% close file if necessary
if beamOutput ~= 1
	fclose(beamOutput);
end

% display cumulative beam angle as additional information
fprintf('Total beam angle = %f sr = %f deg^2.\n', totalSteradians, radtodeg(radtodeg(totalSteradians)));

end


function angleDegree = radtodeg( angleRadiant )

angleDegree = angleRadiant / pi * 180;
end
