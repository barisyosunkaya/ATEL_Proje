%% build_synthetic_fhss_stream.m
% MPACT tekil burst dosyalarini sentetik frekans kaydirma ve cerceve araligiyla
% birlestirerek Pluto SDR uyumlu (cf32) surekli bir FHSS IQ dosyasi uretir.

clear; clc; close all;

%% 1. Ayarlar
dataDir    = '/Users/zeynepsahin/Desktop/feature_ext/Futaba_T8FG'; % MAT dosyalarinin yolu
outputIq   = fullfile(pwd, 'synthetic_futaba_fhss_stream.iq');     % Pluto SDR icin cikti dosyasi

Fs_sdr     = 20e6;       % Pluto SDR ornekleme hizi (20 MSa/s)
bwDDC      = 20e6;       % Filtre bant genisligi (20 MHz)
T_frame    = 4e-3;       % Kumanda cerceve periyodu (4 ms = 0.25 ms burst + 3.75 ms sessizlik)
N_frame    = round(T_frame * Fs_sdr); % 1 cercevedeki toplam ornek sayisi (80,000 ornek)

numHops    = 20;         % Birlestirilecek toplam hop / dosya sayisi

% 20 MHz taban bant (-10 MHz ile +10 MHz) icinde sentetik hop frekanslari (Hz):
% Kenar filtrelerine carpmamasi icin -8 MHz ile +8 MHz arasinda secilir:
hopTable = [-6e6, +2e6, -3e6, +5e6, -1e6, +7e6, -4e6, +3e6, -7e6, +1e6, ...
            -5e6, +4e6, -2e6, +6e6,  0e6, +8e6, -8e6, +2e6, -4e6, +5e6];

%% 2. Dosyalari Listele
files = dir(fullfile(dataDir, '*.mat'));
assert(~isempty(files), 'Klasorde MAT dosyasi bulunamadi!');
[~, ord] = sort({files.name}); files = files(ord);
numHops = min(numHops, numel(files));

fprintf('Sentetik akis olusturuluyor: %d hop, Cerceve suresi = %.1f ms, Toplam sure = %.1f ms\n', ...
    numHops, T_frame*1e3, numHops*T_frame*1e3);

%% 3. Dosyalari Isle, Frekans Kaydir ve Birlestir
fullStream = complex(zeros(numHops * N_frame, 1));

for k = 1:numHops
    fn = fullfile(files(k).folder, files(k).name);
    fprintf('[%2d/%2d] Isleniyor: %s -> Hedef Frekans: %+6.2f MHz\n', ...
        k, numHops, files(k).name, hopTable(k)/1e6);
    
    % a) Ham veriden tasiyiciyi bul
    m = matfile(fn); d = m.data;
    x_raw = double(d.Data(round(0.57*end):end)) * double(d.YInc);
    [P, f] = pwelch(x_raw - mean(x_raw), hann(2^16), 2^15, 2^16, 1/double(d.XInc));
    band = (f >= 2.39e9) & (f <= 2.49e9);
    fb = f(band); Pb = P(band);
    Ps = movmean(Pb, 5); med = median(Ps);
    [~, ip] = max(Ps);
    near = abs(fb - fb(ip)) <= 15e6;
    w = max(Pb(near) - 1.5*med, 0);
    fcEst = sum(fb(near) .* w) / sum(w);
    
    % b) 0 Hz taban banda cevir (mpact_emulate_sdr_iq kullanarak)
    info = mpact_emulate_sdr_iq(fn, '', fcEst, Fs_sdr, bwDDC, ...
        'useVolts', true, 'normalize', false, 'returnIQ', true, 'verbose', false, ...
        'stopAtten', 70, 'lpfRollOff', 0.9);
    y_bb = info.iq(:);
    
    % c) Yalnizca aktif burst kesitini al (~0.25 ms / ~5000 ornek)
    % Filtre gecikmesi payi birakilarak ilk 5000 ornek alinir:
    nBurst = min(5000, numel(y_bb));
    burst  = y_bb(1:nBurst);
    
    % d) SENTETIK MIKSER: Burst'u hop frekansina kaydir
    t = (0:nBurst-1).' / Fs_sdr;
    delta_f = hopTable(mod(k-1, numel(hopTable)) + 1);
    burst_shifted = burst .* exp(1j * 2 * pi * delta_f * t);
    
    % e) Cerceveye yerlestir (0.25 ms sinyal + geriye kalan kisim sessizlik/gurultu)
    frame = complex(zeros(N_frame, 1));
    % İstenirse aradaki sessizlige cok hafif taban gurultusu eklenebilir:
    frame = frame + (randn(N_frame,1) + 1j*randn(N_frame,1)) * 1e-4;
    frame(1:nBurst) = frame(1:nBurst) + burst_shifted;
    
    % f) Ana akisa ekle
    idxStart = (k-1)*N_frame + 1;
    idxEnd   = k*N_frame;
    fullStream(idxStart:idxEnd) = frame;
end

%% 4. Global Normalizasyon ve Diske Kayit (Pluto SDR cf32 formati)
peak = max(abs(fullStream));
if peak > 0
    fullStream = fullStream * (0.8 / peak);
end

fid = fopen(outputIq, 'wb');
assert(fid > 0, 'Cikti dosyasi olusturulamadi!');
tmp = [single(real(fullStream)) single(imag(fullStream))].';
fwrite(fid, tmp(:), 'float32');
fclose(fid);
fprintf('Tamamlandi! Pluto SDR icin dosya hazir: %s\n', outputIq);

%% 5. Birlestirilmis Akisin Spektrogramini Cizdir (FHSS Waterfall Dogrulamasi)
figure('Color', 'w', 'Position', [100 100 1200 600]);
nfft = 1024;
spectrogram(fullStream, hann(nfft), round(0.75*nfft), nfft, Fs_sdr, 'centered', 'yaxis');
colormap('turbo');
title(sprintf('Birlestirilmis Sentetik FHSS Akisi (Toplam Sure: %.1f ms - %d Hop)', (numel(fullStream)/Fs_sdr)*1e3, numHops));
xlabel('Zaman (ms)');
ylabel('Frekans (MHz)');