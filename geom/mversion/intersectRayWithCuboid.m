function [hit, tn, tf] = intersectRayWithCuboid(sourcePos, rayDir, cubMin, cubMax)
% [hit, tn, tf] = intersectRayWithCuboid(sourcePos, rayDir, cubMin, cubMax)
%
% Copyright (c) 2010 by Andreas Keil, TU Munich.

% this code should work in any number of dimensions; get the dimension here
dims = length(sourcePos);

% init near and far intersection distance
tn = -Inf;
tf = Inf;

% intersect with bounding box' x, y (and z) planes by looping over them
for i = 1:dims
	
	% check whether ray is parallel to i planes
	if abs(rayDir(i)) < eps
		
		% check whether ray is outside i planes
		if sourcePos(i) < cubMin(i) || sourcePos(i) > cubMax(i)
			hit = false;
			return;
		end
		
	else
		
		% compute the two intersections
		t1 = (cubMin(i)-sourcePos(i)) / rayDir(i);
		t2 = (cubMax(i)-sourcePos(i)) / rayDir(i);
		
		% sort t1 and t2 so that t1 is intersection with near plane, t2 with far plane
		if (t1 > t2)
			% swap the two values if necessary
			ttmp = t1;
			t1 = t2;
			t2 = ttmp;
		end
		
		% keep largest near intersection
		if t1 > tn, tn = t1; end
		
		% keep smallest far intersection
		if t2 < tf, tf = t2; end
		
		% check if box was missed
		if tn > tf
			hit = false;
			return;
		end
		
		% check if box is behind the ray
		if tf < 0
			hit = false;
			return;
		end
		
	end % if ray is not parallel to the two i planes
	
end % for all plane orientations

% survived all tests, bounding box was hit by ray
hit = true;
