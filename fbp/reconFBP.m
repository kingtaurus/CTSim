function img = reconFBP( sino, geom, window, viewUpsamplingRate, crop )
% Filtered back-projection of CT recontruction
%   This function will call different FBP method based on the geometry
%   parameter of the data for circular scan.
%
%   Now the methods include: fan beam and cone beam circular CT
%                            2D and 3D both
%                            flat panel and curve detector
%                            circular scan only
%
%   Helical scan should call function reconHelical(), because there are
%   many different vesions and more parameters need to be decided.
%
%   For helical CT, 3 methods have been implemented, but this function only call
%   one of it. More details see implementations in the /fbp/ folder.
%
%   There are multiple window functions for ramp filtering can be applied.
%   More details about the window function, see
%           function filt = designFilter2(gamma, window, len, crop)
%           function filt = designEquiAgularFilter2(gamma, window, len, crop)
%
% input:
%       sino    - log sinogram
%       geom    - geometry parameters
%       window  - window function ( default 'hamming', other options: 'ram-lak', 'shepp-logan', 'cosine'  )
%       viewUpsamplingRate (default 1)
%       crop    - frequency crop ratio (default 1)
% output:
%       img     - reconstructed attenuation image (1/cm)
%
% Meng Wu, Stanford University, 2013-04
% Modified 2014-03, the water correction is moved out of the FBP
% reconstruction
% Modified 2014-04, the helical reconstruction is move to reconHelical()

if nargin < 3
    window = 'hamming';
end

if nargin < 4
    viewUpsamplingRate = 1;
end

if nargin < 5
    crop = 1;
end



% view upsampling
if viewUpsamplingRate > 1
    [sino, geom ] = upsamplingViews( sino, geom, viewUpsamplingRate );
elseif viewUpsamplingRate < 0
    %Fourier based downsampling with IIR filter for anti-aliasing
    % [sino, geom] = binningDetectorColumns( sino, geom, -viewUpsamplingRate );
    % Rectangular funcation anti-aliasing and simple downsampling (much faster)
    [sino, geom] = binningDetectorPixels( sino, geom, [-viewUpsamplingRate, 1] );
end

if geom.helicalScan
    
    fprintf('The please use function reconHelical() for more reconstruction options. \n');
    img = reconHelical( sino, geom, window );
    
else
    if length( geom.reconSize ) == 3 && length( geom.detSize ) == 2
        
        geom = editNoViews( geom, size( sino, 3) );
        
        sino = binningDetecotRows( sino, round(  geom.reconSpacing(3) / geom.detSpacing(2) ) );
        img = reconFDK( sino, geom, window, crop );
    elseif length( geom.reconSize ) == 2 && length( geom.detSize ) == 1
        img = reconFBP2d( sino, geom, window, crop );
    end
end

%end of file

end



function geom = editNoViews( geom, noViews )

if geom.noViews == noViews
    return;
end

geom.noViews = noViews;

geom.betas     = geom.betas( 1:noViews );
geom.couchZ    = geom.couchZ(  1:noViews  );

if isfield( geom, 'SADPerProjection' )
    geom.SADPerProjection = geom.SADPerProjection(  1:noViews  );
end

if isfield( geom, 'SDDPerProjection' )
    geom.SDDPerProjection = geom.SDDPerProjection(  1:noViews  );
end

if isfield( geom, 'ADDPerProjection' )
    geom.ADDPerProjection = geom.ADDPerProjection(  1:noViews  );
end

if isfield( geom, 'detOffsetsPerProjection' )
    geom.detOffsetsPerProjection = geom.detOffsetsPerProjection( :,  1:noViews  );
end

if isfield( geom, 'cameraPositions' )
    geom.cameraPositions = geom.cameraPositions( :,  1:noViews  );
end

if isfield( geom, 'PMat' )
    geom.PMat = geom.PMat( :, :,  1:noViews  );
end


end

