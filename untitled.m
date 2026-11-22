filename = "C:\Users\User\Downloads\YUNZHUO H16\YUNZHUO H16\pack1_1-2s.iq";

fs = 100e6; % 100 MSPS sampling rate

% 1. Load the entire IQ file in single precision (800 MB for 100M samples)
fid = fopen(filename, 'rb');
assert(fid ~= -1, 'Could not open file.');
raw_data = fread(fid, Inf, 'float32=>single');
fclose(fid);

iq_signal = raw_data(1:2:end) + 1j * raw_data(2:2:end);
clear raw_data;

num_samples = length(iq_signal);
duration_sec = num_samples / fs;
fprintf('Loaded %d complex samples (%.3f s total)\n', num_samples, duration_sec);

figure('Name', 'Complete 1.00 s IQ Signal Overview', 'Position', [80, 60, 1200, 820]);

% =========================================================================
% SUBPLOT 1: Full 1.00 s I & Q Waveforms (Peak-Preserving Min/Max Binning)
% =========================================================================
num_bins = 20000; % 50 us per bin across 1.00 s (preserves 100% of peaks)
samples_per_bin = floor(num_samples / num_bins);
usable_samples = num_bins * samples_per_bin;

I_mat = reshape(real(iq_signal(1:usable_samples)), samples_per_bin, num_bins);
Q_mat = reshape(imag(iq_signal(1:usable_samples)), samples_per_bin, num_bins);

% Interleave [min; max] per bin so every spike across all 100M samples is drawn
I_env = reshape([min(I_mat, [], 1); max(I_mat, [], 1)], 1, []);
Q_env = reshape([min(Q_mat, [], 1); max(Q_mat, [], 1)], 1, []);
t_env = repelem(((0:num_bins-1) + 0.5) * (samples_per_bin / fs), 2);
clear I_mat Q_mat;

ax1 = subplot(3, 1, 1);
plot(t_env, I_env, 'b', 'DisplayName', 'In-Phase (I)');
hold on;
plot(t_env, Q_env, 'r', 'DisplayName', 'Quadrature (Q)');
hold off;
title(sprintf('Full %.2f s Time-Domain I & Q Waveforms (All %dM Samples, Peak-Preserved)', ...
    duration_sec, round(num_samples/1e6)));
xlabel('Time (s)');
ylabel('Amplitude');
legend('Location', 'northeast');
grid on;
xlim([0, duration_sec]);

% =========================================================================
% SUBPLOT 2: Full 1.00 s Spectrogram (Time vs. Frequency vs. Power)
% =========================================================================
n_fft = 1024;
num_slices = 2000; % 500 us time step across the full 1.00 s capture
step_samples = floor((num_samples - n_fft) / (num_slices - 1));
win = single(0.5 * (1 - cos(2 * pi * (0:n_fft-1)' / (n_fft - 1)))); % Hann window

% Extract 2000 frames evenly across the entire 1.00 s file and compute STFT
idx_starts = 1 + (0:num_slices-1) * step_samples;
frames = zeros(n_fft, num_slices, 'single');
for k = 1:num_slices
    frames(:, k) = iq_signal(idx_starts(k) : idx_starts(k) + n_fft - 1) .* win;
end

stft_mag = abs(fftshift(fft(frames, n_fft, 1), 1)) / sum(win);
stft_db = 20 * log10(stft_mag + eps);
clear frames;

t_spec = (idx_starts - 1 + n_fft/2) / fs;
f_mhz = ((-n_fft/2 : n_fft/2 - 1) * (fs / n_fft)) / 1e6;

ax2 = subplot(3, 1, 2);
imagesc(t_spec, f_mhz, stft_db);
axis xy;
colormap(jet);
cb = colorbar;
cb.Label.String = 'dBFS';
clim([max(stft_db(:)) - 55, max(stft_db(:))]);
title('Full 1.00 s Spectrogram (All Bursts Across Time & Frequency)');
xlabel('Time (s)');
ylabel('Frequency Offset (MHz)');

% Link time axes so zooming in Subplot 1 automatically zooms Subplot 2
linkaxes([ax1, ax2], 'x');

% =========================================================================
% SUBPLOT 3: Full-Capture Max-Hold & Average Spectrum
% =========================================================================
psd_avg = mean(stft_db, 2);
psd_max = max(stft_db, [], 2);

subplot(3, 1, 3);
plot(f_mhz, psd_max, 'r', 'LineWidth', 0.9, 'DisplayName', 'Max-Hold (Peak Bursts)');
hold on;
plot(f_mhz, psd_avg, 'k', 'LineWidth', 0.9, 'DisplayName', '1.00 s Time-Average');
hold off;
title('Full-Capture Frequency Spectrum (1.00 s Average vs. Peak Max-Hold)');
xlabel('Frequency Offset from LO (MHz)');
ylabel('Magnitude (dBFS)');
legend('Location', 'northeast');
grid on;
xlim([-fs/(2e6), fs/(2e6)]);