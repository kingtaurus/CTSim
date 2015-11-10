function I = paintEllipse(I, center, semiaxes, zrot, val)

rotationMatrix = [cos(zrot), -sin(zrot) 0; sin(zrot), cos(zrot), 0; 0, 0, 1];
extentsXY = max(semiaxes(1:2));
extentsZ = semiaxes(3);
minIter = floor(center-[extentsXY extentsXY extentsZ]);
maxIter = ceil(center+[extentsXY extentsXY extentsZ]);
for x = minIter(1):maxIter(1)
	for y = minIter(2):maxIter(2)
		for z = minIter(3):maxIter(3)
			if isInsideEllipse([x, y, z], center, semiaxes, rotationMatrix)
				I(x, y, z) = val;
			end
		end
	end
end
