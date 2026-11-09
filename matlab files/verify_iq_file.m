%% verify_iq_file_safe.m
% Hatasız IQ doğrulama: pencere/overlap/nfft sinyal uzunluğuna göre otomatik ayarlanır.

clear; clc; close all;

%% ---- Kullanıcı ayarları ----
iqFile = "dji_rc_24g_sc16.iq";   % <-- dosyan
Fs     = 20e6;                  % <-- SDR Fs
fmt    = "sc16";                % "sc16" veya "cf32"
Nread  = 2e6;                   % maksimum okunacak IQ sample (dosya kısa ise otomatik azalır)

%% ---- Dosya boyutundan maksimum okunabilir sample sayısını hesapla ----
d = dir(iqFile);
assert(~isempty(d), "Dosya bulunamadı: %s", iqFile);
bytes = d.bytes;

switch lower(fmt)
    case "sc16"
        bytesPerIQ = 4;   % I(int16)+Q(int16) = 2*2 bytes
    case "cf32"
        bytesPerIQ = 8;   % I(float32)+Q(float32) = 2*4 bytes
    otherwise
        error("fmt sc16 veya cf32 olmalı.");
end

Ntotal = floor(bytes / bytesPerIQ);
Nread  = min(Nread, Ntotal);

fprintf("File bytes: %d | Total IQ: %d | Reading: %d\n", bytes, Ntotal, Nread);

%% ---- IQ oku ----
x = read_iq_interleaved(iqFile, fmt, Nread);

L = length(x);
fprintf("Read IQ length L = %d\n", L);
assert(L >= 256, "Sinyal çok kısa (L=%d). Nread'i artır veya dosya/formatı kontrol et.", L);

%% ---- Basit istatistik ----
dc = mean(x);
fprintf("DC mean = %.3e + j%.3e\n", real(dc), imag(dc));
fprintf("Max|x|=%.3f, RMS|x|=%.3f\n", max(abs(x)), rms(abs(x)));

%% ---- Time domain plot ----
Ns = min(20000, L);
t = (0:Ns-1)/Fs;

figure;
plot(t*1e3, real(x(1:Ns))); grid on;
xlabel("Time (ms)"); ylabel("I");
title("Time domain (I)");

figure;
plot(t*1e3, imag(x(1:Ns))); grid on;
xlabel("Time (ms)"); ylabel("Q");
title("Time domain (Q)");

%% ---- Güvenli PSD (pwelch) ----
% Pencere uzunluğunu L'e göre seç
winL = min(8192, max(256, 2^floor(log2(L/8))));  % ~L/8'e yakın power-of-2, [256..8192]
winL = min(winL, L);                              % asla L'den büyük değil
nover = min(floor(0.5*winL), winL-1);             % overlap < winL garantisi
nfft = min(8192, 2^floor(log2(L)));               % nfft <= L ve <=8192

figure;
[pxx,f] = pwelch(x([1:end/2]), hann(winL), nover, nfft, Fs, "centered");
plot(f/1e6, 10*log10(pxx)); grid on;
xlabel("Frequency (MHz)"); ylabel("PSD (dB/Hz)");
title(sprintf("PSD (Welch) | L=%d, win=%d, ov=%d, nfft=%d", L, winL, nover, nfft));

%% ---- Güvenli Spectrogram ----
winS = min(2048, max(256, 2^floor(log2(L/20))));  % ~L/20
winS = min(winS, L);
ovS  = min(floor(0.75*winS), winS-1);
nfftS = min(4096, 2^floor(log2(max(256, winS)))); % nfftS >= winS değil şart değil ama makul

figure;
spectrogram(x, hann(winS), ovS, nfftS, Fs, "centered", "yaxis");
title(sprintf("Spectrogram | win=%d, ov=%d, nfft=%d", winS, ovS, nfftS));
colorbar;

%% ---- I-Q scatter ----
K = min(200000, L);
figure;
scatter(real(x(1:K)), imag(x(1:K)), 3, '.'); grid on;
axis equal;
xlabel("I"); ylabel("Q");
title("I-Q scatter");

fprintf("OK: doğrulama tamam.\n");

%% ---- Helper: read interleaved IQ ----
function x = read_iq_interleaved(filename, fmt, Nread)
    fid = fopen(filename,'rb');
    assert(fid>0, "Dosya açılamadı: %s", filename);
    c = onCleanup(@() fclose(fid));

    switch lower(fmt)
        case "sc16"
            raw = fread(fid, 2*Nread, 'int16=>double');
            raw = raw(:);
            if numel(raw) < 2
                x = complex([]);
                return;
            end
            raw = raw / 32768;
            I = raw(1:2:end);
            Q = raw(2:2:end);
            x = complex(I, Q);

        case "cf32"
            raw = fread(fid, 2*Nread, 'float32=>double');
            raw = raw(:);
            if numel(raw) < 2
                x = complex([]);
                return;
            end
            I = raw(1:2:end);
            Q = raw(2:2:end);
            x = complex(I, Q);

        otherwise
            error("Bilinmeyen fmt: %s", fmt);
    end
end


