function imgVar = reconFBPVar( sino, geom, window, type, crop )
% Filtered back-projection of CT recontruction's variance
%   This function will call different FBP method based on the geometry
%   parameter of the data for circular scan.
%
% input:
%       sino    - log sinogram
%       geom    - geometry parameters
%       window  - window function ( default 'hamming', other options: 'ram-lak', 'shepp-logan', 'cosine'  )
%       viewUpsamplingRate (default 1)
%       crop    - frequency crop ratio (default 1)
% output:
%       imgVar     - reconstructed attenuation image (1/cm)
%
% Meng Wu, FAU , 2014

if nargin < 3
    window = 'hamming';
end

if nargin < 4
    type = 'all';
end

if nargin < 5
    crop = 1;
end

if geom.helicalScan
    fprintf('Error: not implemented for helical scan. \n');
else
    if length( geom.reconSize ) == 3 && length( geom.detSize ) == 2
        
        nbin = ceil(  geom.reconSpacing(3) / geom.detSpacing(2) );
        sino = binningDetecotRows( sino, nbin ) / nbin ;
        imgVar = reconFDKVar( sino, geom, window, crop, type );
       
    elseif length( geom.reconSize ) == 2 && length( geom.detSize ) == 1
       
        fprintf('Error: not implemented for 2D scan. \n');
        imgVar = reconFBP2d( sino, geom, window, crop );
        
    end
    
end


%imgVar(imgVar < 0) = 0;
%end of file