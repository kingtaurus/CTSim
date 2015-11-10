function img = reconHelical( sino, geom, window, reconMethod, segmentLength, tiltedFilter,  viewUpsamplingRate, crop )
% Filtered back-projection for helical X-ray CT recontruction
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
%       segmentLength 
%       tiltedFilter
%       crop    - frequency crop ratio (default 1)
% output:
%       img     - reconstructed attenuation image (1/cm)
%
% Meng Wu, Stanford University, 2013-04


if nargin < 3,   window = 'hamming';end
if nargin < 4,   reconMethod = 'wfdk'; end
if nargin < 5,   segmentLength = 'short'; end
if nargin < 6,   tiltedFilter = false; end
if nargin < 7,   viewUpsamplingRate = 1; end
if nargin < 8,   crop = 1; end

fprintf('FBP Hilical CT  reconstruction with %s line segement and %s method:  \n', segmentLength, reconMethod);


%% calculate number of projection in each segment
if geom.flatPanel
    fanAngle    = 2 * atan( geom.detSpacing(1) * geom.detSize(1) / ( 2 * geom.SDD ) );
else
    fanAngle    = geom.detSpacing(1) * geom.detSize(1) / geom.SDD ;
end

if strcmpi(segmentLength, '2pi')   % 2 pi segments
    factorOfTurns = 1;
    fprintf('\tUse 2 pi rotation for reconstruction. \n');
else % pi + fan
    factorOfTurns = ( pi + fanAngle + 0.2 ) / ( 2 * pi );
    fprintf('\tUse pi + fan rotation for reconstruction. \n');
end
% detect not sufficient slice coverage caused by the high pitch
if geom.pitch * factorOfTurns > ( geom.SAD - geom.reconSize(1) * geom.reconSpacing(1) / 2 ) / geom.SDD
    fprintf( '\tWarning: the pitch may be too large to coverage entire FOV. \n');
end

noViewsPerSlice =  ceil( factorOfTurns * geom.noViewsTurn / 2 ) * 2;

%% view upsampling
if viewUpsamplingRate > 1
    [sino, geom ] = upsamplingViews( sino, geom, viewUpsamplingRate );
elseif viewUpsamplingRate < 0
    %Fourier based downsampling with IIR filter for anti-aliasing
    % [sino, geom] = binningDetectorColumns( sino, geom, -viewUpsamplingRate );
    % Rectangular funcation anti-aliasing and simple downsampling (much faster)
    [sino, geom] = binningDetectorPixels( sino, geom, [-viewUpsamplingRate, 1] );
end

%% ramp filter
if geom.flatPanel
    H = geom.detSpacing(1) * designFilter2( geom.detSpacing(1), window, geom.detSize(1), crop);
else
    H = geom.detSpacing(1) * designEquiangularFilter2(geom.detSpacing(1), geom.SDD, window, geom.detSize(1), crop);
end

%% actual recon here
if geom.helicalScan
    
    if strcmpi(reconMethod, 'fbp')
        img = reconHelicalFBP( sino, geom, window, noViewsPerSlice, tiltedFilter, crop );
    elseif strcmpi(reconMethod, 'fdk')
        img = reconHelicalFDK( sino, geom, H, noViewsPerSlice, tiltedFilter );
    elseif strcmpi(reconMethod, 'wfdk')
        img = reconHelicalWeightedFDK(  sino, geom, H, noViewsPerSlice, tiltedFilter );
    elseif strcmpi(reconMethod, 'awfdk')
        img = reconHelicalApproxWeightedFDK(  sino, geom, H, noViewsPerSlice, tiltedFilter );
    elseif strcmpi(reconMethod, 'vwfdk')
        img = reconHelicalViewWeightedFDK(  sino, geom, H, noViewsPerSlice, tiltedFilter );
    else
        error( 'Error: unknown helical reconstruction method. \n');
    end
    
else
    fprintf('The please use function reconHelical() for more reconstruction options. \n');
    img = reconFBP( sino, geom, window, 1, crop );
end


img(img < 0) = 0;

fprintf('\nDone.\n');

