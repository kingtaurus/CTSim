function filt = designWhindow2( window, len, crop)

% this kernel design is based on pp74 of Kak & Slaney 's book. It works
% much better than the other window (designFilter.m)

order = max(64,2^nextpow2(2*len-1));

index = -order/2+1:order/2-1;

h = zeros(size(index));
h(index==0) = 1;

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