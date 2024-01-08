function [y2] = image_BandpassFilt(y,hpf,lpf,fs)
%image_BandpassFilt Band pass filtering on image space
%   chunk of code from hmrR_BandpassFilt.m
y2 = y;

% Check for NaN
if max(max(isnan(y)))
    warning('Input to hmrR_BandpassFilt contains NaN values. Add hmrR_PreprocessIntensity_NAN to the processing stream.');
    return
end
% Check for finite values
if max(max(isinf(y)))
    warning('Input to hmrR_BandpassFilt must be finite.');
    return
end

% Check that cutoff < nyquist
if lpf / (fs / 2) > 1 || hpf / (fs / 2) > 1
    warning(['hmrR_BandpassFilt cutoff cannot exceed the folding frequency of the data with sample rate ', num2str(fs), ' hz.']);
    return
end

% low pass filter
lpf_norm = lpf / (fs / 2);
if lpf_norm > 0  % No lowpass if filter is
    FilterOrder = 3;
    [z, p, k] = butter(FilterOrder, lpf_norm, 'low');
    [sos, g] = zp2sos(z, p, k);
    y2 = filtfilt(sos, g, double(y));
end

% high pass filter
hpf_norm = hpf / (fs / 2);
if hpf_norm > 0
    FilterOrder = 5;
    [z, p, k] = butter(FilterOrder, hpf_norm, 'high');
    [sos, g] = zp2sos(z, p, k);
    y2 = filtfilt(sos, g, y2);
end

%data2(ii).SetDataTimeSeries(y2);
end