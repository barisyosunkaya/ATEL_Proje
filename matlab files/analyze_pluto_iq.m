%% analyze_pluto_iq.m
% Analyzes an ADALM-PLUTO SDR received IQ capture (20 MS/s) using the
% existing extract_rf_features.m pipeline and burst/hop analysis logic.
%
% Key adaptations from the 100 MS/s extract_features_batch.m pipeline:
%   1) No second decimation (Fs is already 20 MS/s; decim = 1).
%   2) STFT window scaled from 1024 (@100 MS/s) to 256 (@20 MS/s) to keep
%      ~12.8 us time resolution and 78.1 kHz bin width.
%   3) AD9363 DC offset removal + optional 0 Hz LO-leakage notch during
%      frame detection so residual LO spikes do not create false bursts.
%   4) Adjusted SNR/excess thresholds for over-the-air TX->RX captures and
%      dynamic usableBW adjustment to prevent spectral wrap-around on offset bursts.
%
% Requires: extract_rf_features.m in the same folder or MATLAB path,
%           Signal Processing Toolbox, MATLAB R2019b+.

clear; clc; close all;

%% ---- 1. Input File & SDR Capture Settings ----
S.iqFile     = '';                 % Full path to Pluto IQ file; '' opens a file picker
S.dataFormat = 'cf32';             % 'cf32' (interleaved float32), 'sc16' (interleaved int16), or 'mat'
S.Fs         = 20e6;               % Pluto RX sample rate (20 MS/s after 100 -> 20 MHz decimation)
S.fc         = 2.44e9;             % Pluto RX center frequency (Hz) — set to your tuned LO frequency
S.deviceTag  = 'Pluto_RX_Capture'; % Label used in output folder/plots (overridden if .xml exists)
S.outRoot    = fullfile(pwd, 'pluto_analysis_out');

%% ---- 2. Analysis Settings (Adapted for 20 MS/s Over-the-Air RX) ----
S.sliceSec     = 1.0;    % Process capture in 1.0 s slices (handles shorter or longer files automatically)
S.win          = 256;    % STFT window at 20 MS/s (12.8 us, 78.1 kHz/bin; matches 1024 @ 100 MS/s)
S.blockSec     = 0.05;   % Block length for noise floor estimation (s)
S.thr_dB       = 10;     % Active-frame threshold above noise floor (dB) (lowered from 12 for OTA RX)
S.edge_dB      = 6;      % Carrier band edge threshold (dB)
S.gapFrames    = 2;      % Empty-frame gap tolerance inside a burst
S.carTol       = 1e6;    % Max carrier drift inside a burst (Hz)
S.chanTol      = 1e6;    % Frequency tolerance for grouping into the same channel (Hz)
S.edgeGuard_Hz = 1.0e6;  % Band-edge guard region at +/- Fs/2 (Hz) -> keeps inner 18 MHz
S.hopMinExc_dB = 14;     % Min burst excess (dB) for hop/feature analysis (lowered from 35 for OTA RX)
S.minBurst_us  = 120;    % Minimum burst duration (us)
S.maxBurstFeat = 40;     % Max bursts extracted per 1 s slice
S.usableBW     = 18e6;   % Usable anti-aliased bandwidth at 20 MS/s (Hz)

% ADALM-PLUTO DC / LO leakage mitigation
S.removeDC     = true;   % Subtract complex mean per block/burst
S.dcNotch_Hz   = 100e3;  % Notch +/-100 kHz around 0 Hz during burst detection (set 0 to disable)

% Same dropped fields as extract_features_batch.m & survey_all_classes.m
S.dropFields = {'dutyCycle','nRuns','nChannels','channelSpacing_Hz','hopSpan_Hz', ...
    'dwell_mean_s','dwell_std_s','hri_mean_s','hri_std_s','hopRate_Hz','gap_mean_s','channelFreqs_Hz'};

%% ---- 3. Select IQ File & Optional XML Metadata ----
if isempty(S.iqFile)
    [fn, fp] = uigetfile({'*.iq;*.raw;*.bin;*.dat;*.cf32;*.cs16;*.mat', 'IQ Files (*.iq, *.raw, *.bin, *.cf32, *.cs16, *.mat)'; ...
                          '*.*', 'All Files'}, 'Select ADALM-PLUTO Received IQ File');
    if isequal(fn, 0), disp('Cancelled.'); return; end
    S.iqFile = fullfile(fp, fn);
end
assert(isfile(S.iqFile), 'IQ file not found: %s', S.iqFile);

[fDir, fBase, fExt] = fileparts(S.iqFile);
if strcmpi(fExt, '.cs16') || strcmpi(fExt, '.sc16')
    S.dataFormat = 'sc16';
elseif strcmpi(fExt, '.mat')
    S.dataFormat = 'mat';
end

% Check if a companion XML exists (only override fc/Drone, keep Fs = 20e6 unless user wants XML Fs)
xmlCand = fullfile(fDir, [fBase '.xml']);
if ~isfile(xmlCand)
    xList = dir(fullfile(fDir, '*.xml'));
    if numel(xList) == 1, xmlCand = fullfile(xList.folder, xList.name); end
end
if isfile(xmlCand)
    meta = readMeta(xmlCand);
    if isfinite(meta.CenterFrequency), S.fc = meta.CenterFrequency; end
    if ~isempty(meta.Drone), S.deviceTag = meta.Drone; end
    fprintf('Loaded metadata from %s (fc = %.4f GHz, Drone = %s)\n', xmlCand, S.fc/1e9, S.deviceTag);
end

outDir = fullfile(S.outRoot, fBase);
if ~isfolder(outDir), mkdir(outDir); end

%% ---- 4. Determine Total Samples & Slice Boundaries ----
Fs = S.Fs;
fc = S.fc;

switch lower(S.dataFormat)
    case 'cf32'
        bytesPerSample = 8;
        d = dir(S.iqFile);
        Ntot = floor(d.bytes / bytesPerSample);
        fid = fopen(S.iqFile, 'rb');
        assert(fid > 0, 'Could not open file: %s', S.iqFile);
        cleaner = onCleanup(@() fclose(fid));
        matIQ = [];
    case 'sc16'
        bytesPerSample = 4;
        d = dir(S.iqFile);
        Ntot = floor(d.bytes / bytesPerSample);
        fid = fopen(S.iqFile, 'rb');
        assert(fid > 0, 'Could not open file: %s', S.iqFile);
        cleaner = onCleanup(@() fclose(fid));
        matIQ = [];
    case 'mat'
        fid = -1;
        matIQ = loadMatIQ(S.iqFile);
        Ntot = numel(matIQ);
    otherwise
        error('Unsupported dataFormat: %s (use cf32, sc16, or mat)', S.dataFormat);
end

assert(Ntot >= 4 * S.win, 'IQ file is too short (%d samples).', Ntot);
totalDur_s = Ntot / Fs;
fprintf('File: %s\nFormat: %s | Fs = %.2f MS/s | fc = %.4f GHz | Samples = %d (%.3f s)\n', ...
    S.iqFile, S.dataFormat, Fs/1e6, fc/1e9, Ntot, totalDur_s);

% Design baseband anti-wrap LPF at Fs = 20 MS/s (passband 9 MHz)
[nT, Wn, beta, ftype] = kaiserord([9e6, 9.8e6], [1 0], [0.01 10^(-60/20)], Fs);
nT = nT + mod(nT, 2);
hLP = fir1(nT, Wn, ftype, kaiser(nT + 1, beta), 'noscale');
trimN = ceil(numel(hLP) / 2) + 10;

sliceSamples = round(S.sliceSec * Fs);
nSlices = max(1, ceil(Ntot / sliceSamples));

%% ---- 5. Process Slices ----
All_Ft = table();
All_H  = table();
All_Bu = table();

for sl = 1:nSlices
    n0 = (sl - 1) * sliceSamples;
    N  = min(sliceSamples, Ntot - n0);
    N  = S.win * floor(N / S.win);
    if N < 4 * S.win
        fprintf('[Slice %d/%d] Remaining samples (%d) < 4*win, skipping tail.\n', sl, nSlices, N);
        break;
    end

    tOffset  = n0 / Fs;
    tEnd     = (n0 + N) / Fs;
    sliceTag = sprintf('%.2f-%.2fs', tOffset, tEnd);
    fprintf('\n[Slice %d/%d] Processing %s (%d samples)...\n', sl, nSlices, sliceTag, N);

    % --- A) Time-frequency activity with STFT ---
    blkN = S.win * floor(S.blockSec * Fs / S.win);
    blkN = max(blkN, 4 * S.win);
    nBlk = ceil(N / blkN);

    act = false(0,1);  car = [];  bwF = [];  exF = [];
    noiseDB = [];  Pshow = [];  fShow = [];

    for b = 1:nBlk
        i0 = n0 + (b - 1) * blkN;
        nb = S.win * floor(min(blkN, N - (b - 1) * blkN) / S.win);
        if nb < S.win, break; end

        x = readIQChunk(fid, matIQ, S.dataFormat, i0, nb);
        if S.removeDC, x = x - mean(x); end

        [Sx, f] = spectrogram(x, hann(S.win), 0, S.win, Fs, 'centered');
        Pdb = 10 * log10(abs(Sx).^2 + eps);

        if isempty(noiseDB)
            sP = sort(Pdb, 2);
            noiseDB = sP(:, max(1, round(0.25 * size(sP, 2)))) + 5.4;
            Pshow = Pdb;  fShow = f;
        end

        [a, c, w, e] = frameCarrier(Pdb, noiseDB, f, S);
        act = [act; a]; car = [car; c]; bwF = [bwF; w]; exF = [exF; e]; %#ok<AGROW>
    end
    nFr = numel(act);

    % --- B) Bursts, channels, hop groups ---
    Bu = segmentBursts(act, car, bwF, exF, S);
    if height(Bu) == 0
        warning('  -> No bursts detected in %s (try lowering S.thr_dB).', sliceTag);
        continue;
    end

    Bu.Slice     = repmat(string(sliceTag), height(Bu), 1);
    Bu.tStart_s  = tOffset + (Bu.sFrame - 1) * S.win / Fs;
    Bu.dur_s     = Bu.nFrames * S.win / Fs;
    Bu.atStart   = Bu.sFrame == 1;
    Bu.atEnd     = Bu.eFrame == nFr;
    Bu.fcAbs_GHz = (fc + Bu.fcRel_Hz) / 1e9;

    minFrHop = ceil(S.minBurst_us * 1e-6 * Fs / S.win) + 2;
    isLong = Bu.nFrames >= minFrHop & Bu.exc_dB >= S.hopMinExc_dB;
    Bu = Bu(isLong, :);
    nB = height(Bu);

    if nB == 0
        warning('  -> No bursts met duration/SNR criteria in %s (max excess = %.1f dB).', ...
            sliceTag, max(exF, [], 'omitnan'));
        continue;
    end

    [ids, centers, counts] = assignChannels(Bu.fcRel_Hz, S.chanTol);
    Bu.chan = ids;
    gs = find([true; diff(ids) ~= 0]);
    ge = [gs(2:end) - 1; nB];
    Gp = table();
    for k = 1:numel(gs)
        t0_g = Bu.tStart_s(gs(k));
        t1_g = Bu.tStart_s(ge(k)) + Bu.dur_s(ge(k));
        Gp = [Gp; table(k, ids(gs(k)), (fc + centers(ids(gs(k))))/1e9, t0_g, t1_g, t1_g - t0_g, ge(k) - gs(k) + 1, ...
            Bu.atStart(gs(k)), Bu.atEnd(ge(k)), ...
            'VariableNames', {'Group','Chan','Chan_GHz','t0_s','t1_s','Dwell_s','nBursts','atStart','atEnd'})]; %#ok<AGROW>
    end

    H = hopSummary(Bu, Gp, centers, counts, act);
    H_tbl = struct2table(H, 'AsArray', true);
    H_tbl.Slice = string(sliceTag);
    H_tbl = movevars(H_tbl, 'Slice', 'Before', 1);
    All_H = [All_H; H_tbl]; %#ok<AGROW>

    % --- C) Per-burst feature extraction via extract_rf_features.m ---
    minFr = ceil(S.minBurst_us * 1e-6 * Fs / S.win) + 2;
    sel = find(Bu.nFrames >= minFr & ~Bu.atStart & ~Bu.atEnd);
    if isempty(sel)
        % Fallback if capture is short and burst touches start/end boundary
        sel = find(Bu.nFrames >= minFr);
    end
    if numel(sel) > S.maxBurstFeat
        sel = sel(unique(round(linspace(1, numel(sel), S.maxBurstFeat))));
    end

    Ft = table();
    for q = sel(:).'
        i1  = n0 + Bu.sFrame(q) * S.win;
        len = (Bu.nFrames(q) - 2) * S.win;
        if len < 4 * 256 + 4 * trimN
            % Include outer frames if interior is slightly shorter than 4*winSTFT
            i1  = n0 + (Bu.sFrame(q) - 1) * S.win;
            len = Bu.nFrames(q) * S.win;
        end

        x = readIQChunk(fid, matIQ, S.dataFormat, i1, len);
        if S.removeDC, x = x - mean(x); end

        % Shift burst to 0 Hz baseband and filter to suppress wrapped out-of-band noise
        xb = x .* exp(-1j * 2 * pi * Bu.fcRel_Hz(q) / Fs * (0:len-1).');
        y  = filtfilt(hLP, 1, xb);

        % Effective unaliased bandwidth after shifting by fcRel_Hz inside a 20 MHz band
        effBW = max(2e6, min(S.usableBW, Fs - 2 * abs(Bu.fcRel_Hz(q))));

        try
            Fb = extract_rf_features(y, Fs, fc + Bu.fcRel_Hz(q), ...
                'trim', trimN, 'winSTFT', 256, 'usableBW', effBW);
            Fb = rmfield(Fb, S.dropFields);
            Ft = [Ft; [table(string(sliceTag), q, Bu.tStart_s(q), Bu.dur_s(q), Bu.fcAbs_GHz(q), ...
                'VariableNames', {'Slice','Burst','tStart_s','dur_s','fcAbs_GHz'}), struct2table(Fb)]]; %#ok<AGROW>
        catch ME
            warning('  -> %s Burst %d feature error: %s', sliceTag, q, ME.message);
        end
    end

    % --- D) 4-panel diagnostic overview plot ---
    fSlice = figure('Color', 'w', 'Position', [50 50 1400 900], 'Visible', 'off');
    try, theme(fSlice, 'light'); catch, end
    tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile;
    nShow = size(Pshow, 2);
    imagesc((tOffset + (0:nShow-1)*S.win/Fs)*1e3, (fc + fShow)/1e9, Pshow);
    axis xy; colormap('jet');
    pk = max(Pshow(:)); clim([pk - 60, pk]); colorbar;
    xlabel('Time (ms)'); ylabel('Frequency (GHz)');
    title(sprintf('Pluto RX Spectrogram (%s - first block)', sliceTag));

    nexttile;
    tF = tOffset + ((0:nFr-1)' + 0.5) * S.win / Fs;
    scatter(tF(act)*1e3, (fc + car(act))/1e9, 8, exF(act), 'filled'); colorbar; grid on;
    xlabel('Time (ms)'); ylabel('Carrier (GHz)');
    title(sprintf('Per-Frame Carrier Trace (%s)', sliceTag));

    nexttile;
    validDwell = Gp.Dwell_s(~Gp.atStart & ~Gp.atEnd) * 1e3;
    if ~isempty(validDwell)
        histogram(validDwell, 20);
    end
    grid on;
    xlabel('Dwell (ms)'); ylabel('Hop count');
    title('Channel Dwell Time (interior groups)');

    nexttile;
    bar((fc + centers)/1e9, counts); grid on;
    xlabel('Channel (GHz)'); ylabel('Burst count');
    title(sprintf('Channel Usage (%d channels)', numel(centers)));

    try, sgtitle(sprintf('%s | %s (Fs = 20 MS/s)', S.deviceTag, fBase), 'Interpreter', 'none'); catch, end
    saveas(fSlice, fullfile(outDir, sprintf('overview_slice_%d.png', sl)));
    close(fSlice);

    All_Ft = [All_Ft; Ft]; %#ok<AGROW>
    All_Bu = [All_Bu; Bu]; %#ok<AGROW>
    fprintf('  -> Done: %d bursts detected, %d bursts analyzed, plot saved: overview_slice_%d.png\n', ...
        nB, height(Ft), sl);
end

%% ---- 6. Aggregate & Print Extracted Features ----
assert(~isempty(All_Ft), 'No valid burst features could be extracted! Check S.fc, S.dataFormat, or lower S.thr_dB / S.hopMinExc_dB.');

featCols = All_Ft.Properties.VariableNames(6:end);

Overall_Feat_Mean = varfun(@(v) mean(v, 'omitnan'), All_Ft(:, 6:end));
Overall_Feat_Mean.Properties.VariableNames = featCols;
Overall_Feat_Mean.TotalSlices = numel(unique(All_Ft.Slice));
Overall_Feat_Mean.TotalBursts = height(All_Ft);
Overall_Feat_Mean = movevars(Overall_Feat_Mean, {'TotalSlices', 'TotalBursts'}, 'Before', 1);

Overall_Feat_Median = varfun(@(v) median(v, 'omitnan'), All_Ft(:, 6:end));
Overall_Feat_Median.Properties.VariableNames = featCols;
Overall_Feat_Median.TotalSlices = numel(unique(All_Ft.Slice));
Overall_Feat_Median.TotalBursts = height(All_Ft);
Overall_Feat_Median = movevars(Overall_Feat_Median, {'TotalSlices', 'TotalBursts'}, 'Before', 1);

hopCols = All_H.Properties.VariableNames(2:end);
Overall_Hop_Mean = varfun(@(v) mean(v, 'omitnan'), All_H(:, 2:end));
Overall_Hop_Mean.Properties.VariableNames = hopCols;
Overall_Hop_Mean.TotalSlices = height(All_H);
Overall_Hop_Mean = movevars(Overall_Hop_Mean, 'TotalSlices', 'Before', 1);

Final_Combined_Vector = [Overall_Hop_Mean, Overall_Feat_Mean(:, 3:end)];

fprintf('\n=================================================================\n');
fprintf('ADALM-PLUTO RX IQ ANALYSIS SUMMARY (%d Slice(s), %d Burst(s))\n', ...
    height(All_H), height(All_Ft));
fprintf('=================================================================\n');

fprintf('\n--- 1. HOP & BURST TIMING SUMMARY (MEAN) ---\n');
disp(rows2vars(Overall_Hop_Mean));

fprintf('\n--- 2. RF & MODULATION FEATURES (MEAN & MEDIAN ACROSS BURSTS) ---\n');
SummaryTbl = table(string(featCols(:)), ...
    table2array(Overall_Feat_Mean(1, 3:end)).', ...
    table2array(Overall_Feat_Median(1, 3:end)).', ...
    'VariableNames', {'Feature', 'Mean', 'Median'});
disp(SummaryTbl);

% Save all CSVs matching extract_features_batch.m output structure
writetable(All_Bu,                fullfile(outDir, 'all_slices_bursts.csv'));
writetable(All_H,                 fullfile(outDir, 'all_slices_hop_summary.csv'));
writetable(All_Ft,                fullfile(outDir, 'all_slices_burst_features.csv'));
writetable(Overall_Feat_Mean,     fullfile(outDir, 'overall_burst_features_mean.csv'));
writetable(Overall_Feat_Median,   fullfile(outDir, 'overall_burst_features_median.csv'));
writetable(Overall_Hop_Mean,      fullfile(outDir, 'overall_hop_summary_mean.csv'));
writetable(Final_Combined_Vector, fullfile(outDir, 'FINAL_COMBINED_FEATURE_VECTOR.csv'));

fprintf('All CSVs and overview figures saved to:\n  %s\n', outDir);

%% ------------------------- Local Helper Functions -------------------------
function x = readIQChunk(fid, matIQ, fmt, startSample, n)
    switch lower(fmt)
        case 'cf32'
            fseek(fid, 8 * startSample, 'bof');
            raw = fread(fid, [2, n], 'float32=>double');
            x = complex(raw(1, :).', raw(2, :).');
        case 'sc16'
            fseek(fid, 4 * startSample, 'bof');
            raw = fread(fid, [2, n], 'int16=>double') / 32768.0;
            x = complex(raw(1, :).', raw(2, :).');
        case 'mat'
            x = matIQ(startSample + 1 : startSample + n);
    end
end

function x = loadMatIQ(matFile)
    L = load(matFile);
    fn = fieldnames(L);
    for k = 1:numel(fn)
        v = L.(fn{k});
        if isnumeric(v) && ~isscalar(v)
            x = double(v(:));
            if isreal(x) && mod(numel(x), 2) == 0
                warning('MAT variable "%s" is real; assuming interleaved I/Q.', fn{k});
                x = complex(x(1:2:end), x(2:2:end));
            end
            return;
        elseif isstruct(v) && isfield(v, 'iq')
            x = double(v.iq(:));
            return;
        end
    end
    error('No IQ vector found inside %s.', matFile);
end

function meta = readMeta(xmlFile)
    txt = fileread(xmlFile);
    numTags = {'SampleRate','CenterFrequency','IFBandwidth','SampleCount','ScaleFactor','ReferenceSNRLevel'};
    for i = 1:numel(numTags)
        tk = regexp(txt, ['<' numTags{i} '>\s*([^<]+?)\s*</' numTags{i} '>'], 'tokens', 'once');
        if isempty(tk), meta.(numTags{i}) = NaN; else, meta.(numTags{i}) = str2double(tk{1}); end
    end
    strTags = {'Drone','DataType','DeviceType'};
    for i = 1:numel(strTags)
        tk = regexp(txt, ['<' strTags{i} '>\s*([^<]+?)\s*</' strTags{i} '>'], 'tokens', 'once');
        if isempty(tk), meta.(strTags{i}) = ''; else, meta.(strTags{i}) = tk{1}; end
    end
end

function [act, car, bw, exc] = frameCarrier(Pdb, noiseDB, f, S)
    E = Pdb - noiseDB;
    if S.dcNotch_Hz > 0
        E(abs(f) <= S.dcNotch_Hz, :) = -Inf;
    end
    [pk, ip] = max(E, [], 1);
    nFr = size(E, 2);  nf = size(E, 1);  df = f(2) - f(1);
    act = pk(:) > S.thr_dB;
    car = nan(nFr, 1);  bw = nan(nFr, 1);  exc = pk(:);
    for k = find(act).'
        lo = ip(k);  while lo > 1  && E(lo-1, k) > S.edge_dB, lo = lo - 1; end
        hi = ip(k);  while hi < nf && E(hi+1, k) > S.edge_dB, hi = hi + 1; end
        w = 10.^(E(lo:hi, k)/10) - 1;  w(w < 0) = 0;
        if sum(w) > 0, car(k) = sum(f(lo:hi) .* w) / sum(w); else, car(k) = f(ip(k)); end
        bw(k) = (hi - lo + 1) * df;
        if abs(car(k)) > max(f) - S.edgeGuard_Hz
            act(k) = false;  car(k) = NaN;  bw(k) = NaN;
        end
    end
end

function Bu = segmentBursts(act, car, bwF, exF, S)
    nFr = numel(act);
    sF = [];  eF = [];  cf = [];  bb = [];  ee = [];
    i = 1;
    while i <= nFr
        if ~act(i), i = i + 1; continue; end
        j = i;
        while j < nFr
            nxt = find(act(j+1 : min(nFr, j + 1 + S.gapFrames)), 1, 'first');
            if isempty(nxt), break; end
            jn = j + nxt;
            if abs(car(jn) - car(j)) > S.carTol, break; end
            j = jn;
        end
        sF(end+1, 1) = i;  eF(end+1, 1) = j; %#ok<AGROW>
        cf(end+1, 1) = median(car(i:j), 'omitnan');   bb(end+1, 1) = median(bwF(i:j), 'omitnan'); %#ok<AGROW>
        ee(end+1, 1) = median(exF(i:j), 'omitnan'); %#ok<AGROW>
        i = j + 1;
    end
    Bu = table(sF, eF, eF - sF + 1, cf, bb, ee, ...
        'VariableNames', {'sFrame','eFrame','nFrames','fcRel_Hz','bw_Hz','exc_dB'});
end

function [ids, centers, counts] = assignChannels(fcs, tol)
    [sv, ord] = sort(fcs(:));
    cl = cumsum([true; diff(sv) > tol]);
    ids = zeros(size(fcs(:)));  ids(ord) = cl;
    centers = accumarray(cl, sv, [], @mean);
    counts  = accumarray(cl, 1);
end

function H = hopSummary(Bu, Gp, centers, counts, act)
    H = struct();
    H.nBursts   = height(Bu);
    H.nGroups   = height(Gp);
    H.nChannels = numel(centers);
    H.dutyCycle = mean(act);
    sc = sort(centers);
    if numel(sc) >= 2
        dsc = diff(sc);  H.channelSpacing_Hz = median(dsc(dsc <= 1.5*min(dsc)));
    else, H.channelSpacing_Hz = NaN; end
    H.hopSpan_Hz = max(centers) - min(centers);
    p = counts / sum(counts);
    if numel(p) > 1, H.chanUsageEntropy = -sum(p .* log2(p)) / log2(numel(p)); else, H.chanUsageEntropy = NaN; end
    ok = ~Gp.atStart & ~Gp.atEnd;
    if any(ok)
        H.dwell_mean_s = mean(Gp.Dwell_s(ok));  H.dwell_std_s = std(Gp.Dwell_s(ok));
        H.dwell_median_s = median(Gp.Dwell_s(ok));
    else
        H.dwell_mean_s = NaN;  H.dwell_std_s = NaN;  H.dwell_median_s = NaN;
    end
    t0 = Gp.t0_s;
    if Gp.atStart(1), t0 = t0(2:end); end
    hri = diff(t0);
    if numel(hri) >= 2
        H.hri_mean_s = mean(hri);  H.hri_std_s = std(hri);  H.hri_median_s = median(hri);
        H.hopRate_Hz = 1 / H.hri_median_s;
    else
        H.hri_mean_s = NaN;  H.hri_std_s = NaN;  H.hri_median_s = NaN;  H.hopRate_Hz = NaN;
    end
    if height(Gp) >= 2
        gap = Gp.t0_s(2:end) - Gp.t1_s(1:end-1);
        H.gap_mean_s = mean(gap);
    else
        H.gap_mean_s = NaN;
    end
    okb = ~Bu.atStart & ~Bu.atEnd;
    if any(okb), H.pktDur_median_s = median(Bu.dur_s(okb)); else, H.pktDur_median_s = NaN; end
    sameGrp = Bu.chan(2:end) == Bu.chan(1:end-1);
    if any(sameGrp)
        dt = diff(Bu.tStart_s);  H.pktInterval_median_s = median(dt(sameGrp));
    else
        H.pktInterval_median_s = NaN;
    end
end