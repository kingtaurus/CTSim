function filt = designEquiangularFilter2(gamma, SDD, window, len, crop)
%function filt = designFilter2(gamma, window, len, crop)
%   Design ramp filter for filter back projection reconstruction
%   Curved detecor case (equiangluar).
%
%   Input:
%       gamma   - angle of detector pixel 
%       window  - window function names
%       len     - length of the detector, zero padding not include
%       crop    - frequency crop ratio
%   Output:
%       filt    - ramp filter in frequency domain, zero padding included.
%   
%
%   Based on Prof. Jeffery Fessler's and Prof. Lei Zhu's implementation   
%   Meng Wu, 2013-2014


order = max(64,2^nextpow2(2*len-1));
index = -order/2+1:order/2-1;

h = zeros(size(index));
h(index==0) = 1/(4*gamma^2);

% note: abs(index)<len, instead of abs(index)<len/2, if needed in 'temp';

temp = find(mod(index,2)==1);
h(temp) = -1./( SDD*pi*sin(index(temp)*gamma/ SDD)).^2;

h(order) = 0; % zero-padding

filt = abs(fft(h));
filt = filt(1:order/2+1);

w = 2*pi*(0:size(filt,2)-1)/order;   % frequency axis up to Nyquist

switch window
    case 'ram-lak'
        % Do nothing
    case 'shepp-logan'
        % be careful not to divide by 0:
        filt(2:end) = filt(2:end) .* (sin(w(2:end)/(2*crop))./(w(2:end)/(2*crop)));
    case 'cosine'
        filt(2:end) = filt(2:end) .* cos(w(2:end)/(2*crop));
    case 'hamming'
        filt(2:end) = filt(2:end) .* (.54 + .46 * cos(w(2:end)/crop));
    case 'hann'
        filt(2:end) = filt(2:end) .*(1+cos(w(2:end)./crop)) / 2;
    case 'hann50'
        crop = 0.5;
        filt(2:end) = filt(2:end) .*(1+cos(w(2:end)./crop)) / 2;
    case 'hann75'
        crop = 0.75;
        filt(2:end) = filt(2:end) .*(1+cos(w(2:end)./crop)) / 2;
    case 'hann80'
        crop = 0.80;
        filt(2:end) = filt(2:end) .*(1+cos(w(2:end)./crop)) / 2;
    case 'B50'
        crop = 0.5;
    case 'B75'
        crop = 0.75;
    case 'B80'
        crop = 0.80;
    case 'C50'
        crop = 0.5;
        filt(2:end) = filt(2:end) .* cos(w(2:end)/(2*crop));
    case 'C75'
        crop = 0.75;
        filt(2:end) = filt(2:end) .* cos(w(2:end)/(2*crop));
    case 'C80'
        crop = 0.80;
        filt(2:end) = filt(2:end) .* cos(w(2:end)/(2*crop));
    otherwise
        error('Invalid window selected.');
end

filt(w>pi*crop) = 0;                      % Crop the frequency response
filt = [filt' ; filt(end-1:-1:2)'];    % Symmetry of the window