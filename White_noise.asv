%RUN_AWGN_ON_IQ  Load a cf32 .iq file, add AWGN at chosen SNRs, plot, save.
clear; close all; clc;

%% ----------Configure ----------
fname    = 'C:\Users\User\Desktop\ATEL(1)\DJI_Phantom4Pro_1_0001_cf32.iq';
Fs       = 20e6;            % sampling rate of the recording (Hz)  <-- set yours
B        = [];              % signal bandwidth in Hz (1.4e6), or [] for whole-band SNR
Nread    = 2e6;             % number of complex samples to read (Inf = whole file)
snr_list = [-10 -5 0 5 10];  % target SNRs in dB
saveOut  = false;           % true -> write noisy files next to the original


%%  Read float32
d = dir(fname);
fprintf('File: %.1f MB = %d complex samples\n', d.bytes/1e6, floor(d.bytes/8));
fid = fopen(fname, 'r');
if fid < 0, error('Cannot open file: %s', fname); end
raw = fread(fid, 2*Nread, 'float32');
fclose(fid);
x = double(raw(1:2:end)) + 1j*double(raw(2:2:end));
clear raw;
fprintf('Read %d samples, mean power = %.3g\n', numel(x), mean(abs(x).^2));

%%  Add noise at each SNR
N   = min(2^16, numel(x));
win = 0.5 - 0.5*cos(2*pi*(0:N-1)'/N);                  % Hann window
f   = (-N/2:N/2-1)' * Fs/N / 1e6;                       % MHz
psd = @(v) 10*log10(fftshift(abs(fft(v(1:N).*win)).^2)/(sum(win.^2)*Fs) + eps);

figure('Name','Spectrum vs SNR');
plot(f, psd(x), 'k', 'LineWidth', 1.2); hold on;
leg = {'Original'};

for k = 1:numel(snr_list)
    if isempty(B)
        [y, ~, info] = add_awgn_snr(x, snr_list(k), 'Seed', k);
    else
        [y, ~, info] = add_awgn_snr(x, snr_list(k), 'Fs', Fs, 'Bandwidth', B, 'Seed', k);
    end
    fprintf('SNR target %+4d dB -> whole-band SNR measured %.2f dB\n', ...
            snr_list(k), info.SNR_measured_dB);

    plot(f, psd(y));  leg{end+1} = sprintf('SNR = %d dB', snr_list(k)); %#ok<SAGROW>

    if saveOut
        [p, n, ~] = fileparts(fname);
        tmp = zeros(2*numel(y), 1, 'single');
        tmp(1:2:end) = real(y);  tmp(2:2:end) = imag(y);
        fo = fopen(fullfile(p, sprintf('%s_snr%+03d.iq', n, snr_list(k))), 'w');
        fwrite(fo, tmp, 'float32');  fclose(fo);
    end
end
grid on; xlabel('Frequency (MHz)'); ylabel('PSD (dB/Hz)');
title('Recorded I/Q with added AWGN'); legend(leg, 'Location', 'northeast');

%% ----------- LOCAL FUNCTION ------------
function [y, noise, info] = add_awgn_snr(x, snr_dB, varargin)
% Add white Gaussian noise to real/complex x for a target SNR.
% Options: 'SignalPower' ('measured' | scalar), 'Bandwidth' (Hz),
%          'Fs' (Hz), 'Seed' (integer).
    p = inputParser;
    addParameter(p, 'SignalPower', 'measured');
    addParameter(p, 'Bandwidth', []);
    addParameter(p, 'Fs', []);
    addParameter(p, 'Seed', []);
    parse(p, varargin{:});
    o = p.Results;

    if ~isempty(o.Bandwidth) && isempty(o.Fs)
        error('''Fs'' is required when ''Bandwidth'' is given.');
    end

    isRow = isrow(x);
    if isRow, x = x(:); end

    if ischar(o.SignalPower)
        Ps = mean(abs(x).^2, 1);
    else
        Ps = o.SignalPower * ones(1, size(x, 2));
    end

    snr_lin = 10^(snr_dB/10);
    if isempty(o.Bandwidth)
        Pn = Ps / snr_lin;
    else
        Pn = Ps / snr_lin * (o.Fs / o.Bandwidth);   % in-band SNR
    end

    if ~isempty(o.Seed), rng(o.Seed); end
    if isreal(x)
        noise = sqrt(Pn) .* randn(size(x));
    else
        noise = sqrt(Pn/2) .* (randn(size(x)) + 1j*randn(size(x)));
    end
    y = x + noise;

    info.SignalPower     = Ps;
    info.NoisePower      = Pn;
    info.SNR_measured_dB = 10*log10(mean(abs(x).^2,1) ./ mean(abs(noise).^2,1));
    if isRow, y = y.'; noise = noise.'; end
end