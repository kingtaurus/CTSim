function[grid3Dout] = mat2Grid3D(inputVolume)
import edu.stanford.rsl.conrad.data.*

if isreal(inputVolume)
    grid3Dout = Grid3D(size(inputVolume,1),size(inputVolume,2),size(inputVolume,3),false);
else
    grid3Dout = Grid3DComplex(size(inputVolume,1),size(inputVolume,2),size(inputVolume,3),false);
end
for t=1:size(inputVolume,3)
    slice = inputVolume(:,:,t);
    if isreal(inputVolume)
        g2d = Grid2D(slice(:),size(inputVolume,1),size(inputVolume,2));
    else
        arr = zeros(length(slice(:))*2,1);
        arr(mod(1:length(arr),2)==1) = imag(slice(:));
        arr(mod(1:length(arr),2)==0) = real(slice(:));
        g2d = Grid2DComplexNEW(arr,size(inputVolume,1),size(inputVolume,2));
    end
    
    grid3Dout.setSubGrid(t-1,g2d);
end
