%% extract_features_batch.m
% Batch version: for every pack (pack1, pack2, ...) found under the selected folder, all 1-second
% slices (pack1_0-1s.iq -> last slice) are found, sorted, processed (4-panel overview per slice,
% burst/hop features) and averaged over the whole recording.
%
%  - Select a folder: a drone folder ("DJI FPV COMBO", "JR PROPO XG14", ...), a sub-folder
%    ("VTSBW=10") or the whole dataset root. All *.iq files below it are found recursively and
%    matched with their <pack>.xml in the same folder.
%  - Every pack gets its own output folder that mirrors the dataset layout and is named after the input:
%        features_out\<drone>\<sub-folder>\<pack>\
%    so results of different drones / settings / packs are never overwritten by each other.
%  - Merged tables for the selected folder (one row per pack) are written to
%        features_out\<selected folder>\_batch_overview.csv and _batch_final_vectors.csv
%    (columns SubPath / Pack identify the source, so different drones can be compared).
%  - The per-pack processing is the same as in the friend's version.
%
% Suggested first runs: S.dryRun = true (only lists what would be processed), then S.iqLimit = 2.
% Requires: extract_rf_features.m (current version), Signal Processing Toolbox, MATLAB R2019b+.
clear; clc; close all;

%% ---- Batch settings ----
S.datasetRoot  = 'C:\Users\ruzga\Documents\Dataset_Bitirme';   % EDIT: parent folder of all drone folders
S.rootDir      = '';       % folder to process; '' opens a folder picker
S.yeniDir      = 'C:\Users\ruzga\OneDrive - Bilkent University\okul\EEE 391 Bitirme\Bitirme_MATLAB';
S.outRoot      = fullfile(S.yeniDir, 'features_out_new');
S.dryRun       = false;    % true: only list the packs/slices that would be processed
S.skipExisting = false;    % true: skip packs that already have a result.mat (resume an interrupted run)
S.iqLimit      = Inf;      % max number of 1 s slices per pack (in time order); 2 for a quick test
S.savePerSlice = false;    % true: also save separate CSVs for every 1 s slice

%% ---- Analysis settings ----
S.tStart       = 0;      % s: start inside each slice file
S.tLen         = 1.0;    % s: length of each slice file
S.win          = 1024;   % STFT window (samples). At 100 MS/s: 10.24 us, bin = 97.7 kHz
S.blockSec     = 0.05;   % block length for memory (s)
S.thr_dB       = 12;     % active-frame threshold (dB)
S.edge_dB      = 6;      % carrier band edge threshold (dB)
S.gapFrames    = 2;      % empty-frame tolerance
S.carTol       = 1e6;    % max carrier change inside a burst (Hz)
S.chanTol      = 1e6;    % frequency tolerance for the same channel (Hz)
S.edgeGuard_Hz = 1.5e6;  % band-edge guard region (Hz)
S.hopMinExc_dB = 35;     % threshold for dominant bursts entering the hop analysis (dB)
S.minBurst_us  = 120;    % minimum burst duration (us)
S.maxBurstFeat = 40;     % max bursts for feature extraction in every 1 s slice
S.decim        = 5;      % 100 -> 20 MS/s
S.fp           = 9e6;    % LPF passband edge (Hz)
S.dropFields = {'dutyCycle','nRuns','nChannels','channelSpacing_Hz','hopSpan_Hz', ...
    'dwell_mean_s','dwell_std_s','hri_mean_s','hri_std_s','hopRate_Hz','gap_mean_s','channelFreqs_Hz'};

%% ---- Select folder and find packs ----
if isempty(S.rootDir)
    S.rootDir = uigetdir(S.datasetRoot, 'Select the folder to process (drone folder, sub-folder or dataset root)');
    if isequal(S.rootDir, 0), disp('Cancelled.'); return; end
end
P = findPairs(S.rootDir, S.datasetRoot, S.iqLimit);
if isempty(P) || height(P) == 0
    error('No .iq file found under %s', S.rootDir);
end
[~, ~, g] = unique(P.folder + "|" + P.pack, 'stable');
nP = max(g);
fprintf('Selected: %s\nFound %d pack(s), %d slice file(s):\n', S.rootDir, nP, height(P));
for p = 1:nP
    r = P(g == p, :);
    fprintf('  %-45s %-12s %d slice(s)%s\n', r.relFolder(1), r.pack(1), height(r), ...
        iff(r.xml(1) == "", '   [NO XML - will be skipped]', ''));
end
if S.dryRun, disp('Dry run: nothing processed.'); return; end

relSel = relToRoot(S.rootDir, S.datasetRoot);
if isempty(relSel), relSel = 'ALL'; end
sumDir = fullfile(S.outRoot, relSel);
if ~isfolder(sumDir), mkdir(sumDir); end

%% ---- Process every pack ----
R = repmat(emptyRes(), nP, 1);
packRel = strings(nP, 1);  packName = strings(nP, 1);
tAll = tic;
for p = 1:nP
    rows = P(g == p, :);
    packRel(p) = rows.relFolder(1);  packName(p) = rows.pack(1);
    outDir  = fullfile(S.outRoot, char(packRel(p)), char(packName(p)));
    resFile = fullfile(outDir, 'result.mat');
    fprintf('\n=================================================================\n');
    fprintf('[%d/%d] %s\\%s  (%d slices)\n', p, nP, packRel(p), packName(p), height(rows));
    if S.skipExisting && isfile(resFile)
        L = load(resFile, 'res');  R(p) = L.res;
        fprintf('  skipped (result.mat exists)\n');
        continue;
    end
    if ~isfolder(outDir), mkdir(outDir); end
    tPack = tic;
    if rows.xml(1) == ""
        res = emptyRes();  res.status = 'no xml';
        res.message = 'No matching .xml file in the folder';
    else
        try
            res = processPack(S, rows, outDir);
        catch ME
            fclose('all');
            res = emptyRes();  res.status = 'error';  res.message = ME.message;
            fprintf(2, '  ERROR: %s\n', ME.message);
        end
    end
    res.elapsed = toc(tPack);
    save(resFile, 'res');
    R(p) = res;
    close all force;
    fprintf('  status: %s (%.0f s)\n', res.status, res.elapsed);
end

%% ---- Merged tables for the selected folder ----
ov = table(packRel, packName, 'VariableNames', {'SubPath','Pack'});
ov.Status        = string({R.status}');
ov.Message       = string({R.message}');
ov.XML_Drone     = string({R.xmlDrone}');
ov.Fs_MHz        = [R.Fs]' / 1e6;
ov.fc_GHz        = [R.fc]' / 1e9;
ov.SlicesFound   = [R.nSlicesFound]';
ov.SlicesUsed    = [R.nSlicesUsed]';
ov.LongBursts    = [R.nBursts]';
ov.FeatBursts    = [R.nFeat]';
ov.SNR_median_dB = [R.snr_med]';
ov.Seconds       = [R.elapsed]';
writetable(ov, fullfile(sumDir, '_batch_overview.csv'));

Vc = cell(nP, 1);
for p = 1:nP
    od = fullfile(S.outRoot, char(packRel(p)), char(packName(p)));
    Vc{p} = tagTable(readIfExists(fullfile(od, 'FINAL_COMBINED_FEATURE_VECTOR.csv')), packRel(p), packName(p));
end
V = stackTables(Vc);
if ~isempty(V), writetable(V, fullfile(sumDir, '_batch_final_vectors.csv')); end

fprintf('\nDone in %.1f min. ok: %d | no xml: %d | error: %d\n', toc(tAll)/60, ...
    sum(strcmp({R.status}, 'ok')), sum(strcmp({R.status}, 'no xml')), sum(strcmp({R.status}, 'error')));
fprintf('Per-pack results: %s\nMerged tables:    %s\n', fullfile(S.outRoot, '<sub-folder>\<pack>'), sumDir);

%% ------------------------- local functions (must stay at the END of the file) -------------------------
function out = iff(c, a, b)
    if c, out = a; else, out = b; end
end

function res = emptyRes()
    res = struct('status', '', 'message', '', 'xmlDrone', '', 'Fs', NaN, 'fc', NaN, ...
        'nSlicesFound', NaN, 'nSlicesUsed', NaN, 'nBursts', NaN, 'nFeat', NaN, 'snr_med', NaN, 'elapsed', NaN);
end

function r = relToRoot(p, root)
    p = regexprep(char(p), '[\\/]+$', '');
    root = regexprep(char(root), '[\\/]+$', '');
    if startsWith(p, root, 'IgnoreCase', true)
        r = char(strip(string(p(numel(root)+1:end)), 'left', filesep));
    else
        [~, r] = fileparts(p);
    end
end

function P = findPairs(rootDir, datasetRoot, iqLimit)
    d = dir(fullfile(rootDir, '**', '*.iq'));
    d = d(~[d.isdir]);
    n = numel(d);
    if n == 0, P = table(); return; end
    folder = strings(n, 1);  name = strings(n, 1);  pack = strings(n, 1);
    t0 = nan(n, 1);  t1 = nan(n, 1);  xml = strings(n, 1);  relFolder = strings(n, 1);
    for i = 1:n
        folder(i) = string(d(i).folder);
        name(i)   = string(d(i).name);
        relFolder(i) = string(relToRoot(d(i).folder, datasetRoot));
        tk = regexp(d(i).name, '^(.*)_(\d+)-(\d+)s\.iq$', 'tokens', 'once', 'ignorecase');
        if isempty(tk)
            [~, pk] = fileparts(d(i).name);  pack(i) = string(pk);  t0(i) = 0;  t1(i) = 1;
        else
            pack(i) = string(tk{1});  t0(i) = str2double(tk{2});  t1(i) = str2double(tk{3});
        end
        cand = fullfile(d(i).folder, [char(pack(i)) '.xml']);
        if isfile(cand)
            xml(i) = string(cand);
        else
            x = dir(fullfile(d(i).folder, '*.xml'));
            if numel(x) == 1, xml(i) = string(fullfile(x.folder, x.name)); end
        end
    end
    P = table(folder, name, pack, t0, t1, xml, relFolder);
    P = sortrows(P, {'folder', 'pack', 't0'});
    if isfinite(iqLimit)
        [~, ~, g] = unique(P.folder + "|" + P.pack, 'stable');
        rk = zeros(n, 1);
        for k = 1:max(g)
            idx = find(g == k);  rk(idx) = 1:numel(idx);
        end
        P = P(rk <= iqLimit, :);
    end
end

function T = readIfExists(f)
    if isfile(f), T = readtable(f); else, T = table(); end
end

function T = tagTable(T, relFolder, pack)
    if isempty(T), return; end
    n = height(T);
    id = table(repmat(relFolder, n, 1), repmat(pack, n, 1), 'VariableNames', {'SubPath', 'Pack'});
    T = [id, T];
end

function T = stackTables(C)
    C = C(~cellfun(@isempty, C));
    if isempty(C), T = table(); return; end
    names = {};
    for i = 1:numel(C), names = [names, C{i}.Properties.VariableNames]; end %#ok<AGROW>
    names = unique(names, 'stable');
    for i = 1:numel(C)
        miss = setdiff(names, C{i}.Properties.VariableNames, 'stable');
        for m = 1:numel(miss), C{i}.(miss{m}) = nan(height(C{i}), 1); end
        C{i} = C{i}(:, names);
    end
    T = vertcat(C{:});
end

function res = processPack(S, rows, outDir)
    res = emptyRes();
    S.outDir = outDir;  S.xmlFile = char(rows.xml(1));
    try
        fj = fopen(fullfile(outDir, 'settings.json'), 'w');
        fprintf(fj, '%s', jsonencode(S));  fclose(fj);
    catch
    end

    %% ---- 1) Metadata ----
    meta = readMeta(S.xmlFile);
    assert(contains(meta.DataType, 'Float', 'IgnoreCase', true), ...
        'Bu script "Complex Float" (cf32) icindir; XML DataType=%s', meta.DataType);
    Fs = meta.SampleRate;
    fc = meta.CenterFrequency;
    assert(isfinite(Fs) && Fs > 0 && isfinite(fc), 'XML is missing SampleRate / CenterFrequency');
    nSl = height(rows);
    res.Fs = Fs;  res.fc = fc;  res.xmlDrone = meta.Drone;  res.nSlicesFound = nSl;

    fprintf('Cihaz: %s | Fs = %.0f MS/s | fc = %.4f GHz\n', meta.Drone, Fs/1e6, fc/1e9);
    fprintf('Bulunan dilim sayisi: %d (%d-%ds -> %d-%ds)\n', ...
        nSl, rows.t0(1), rows.t1(1), rows.t0(end), rows.t1(end));

    % Design the LPF once
    Fs2 = Fs / S.decim;
    [nT, Wn, beta, ftype] = kaiserord([S.fp, Fs2 - S.fp], [1 0], [0.01 10^(-70/20)], Fs);
    nT = nT + mod(nT, 2);
    hLP = fir1(nT, Wn, ftype, kaiser(nT + 1, beta), 'noscale');
    trimN = ceil(numel(hLP) / S.decim) + 20;

    %% ---- 2) Loop over all slices ----
    All_Ft = table();
    All_H  = table();
    All_Bu = table();

    for idx = 1:nSl
        nm       = char(rows.name(idx));
        curFile  = fullfile(char(rows.folder(idx)), nm);
        tOffset  = rows.t0(idx); % for the global time stamp
        sliceTag = sprintf('%d-%ds', rows.t0(idx), rows.t1(idx));

        fprintf('[%d/%d] Isleniyor: %s ...\n', idx, nSl, nm);

        d = dir(curFile);
        Ntot = floor(d.bytes / 8);
        n0 = round(S.tStart * Fs);
        N  = min(round(S.tLen * Fs), Ntot - n0);
        N  = S.win * floor(N / S.win);

        fid = fopen(curFile, 'rb');
        if fid <= 0
            warning('Dosya acilamadi, atlaniyor: %s', nm);
            continue;
        end

        % --- A) Time-frequency activity with STFT ---
        blkN = S.win * floor(S.blockSec * Fs / S.win);
        nBlk = ceil(N / blkN);
        act = false(0,1);  car = [];  bwF = [];  exF = [];
        noiseDB = [];  Pshow = [];  fShow = [];

        for b = 1:nBlk
            i0 = n0 + (b-1)*blkN;
            nb = S.win * floor(min(blkN, N - (b-1)*blkN) / S.win);
            if nb < S.win, break; end
            x = readCF32(fid, i0, nb);
            [Sx, f] = spectrogram(x, hann(S.win), 0, S.win, Fs, 'centered');
            Pdb = 10*log10(abs(Sx).^2 + eps);
            if isempty(noiseDB)
                sP = sort(Pdb, 2);
                noiseDB = sP(:, max(1, round(0.25*size(sP,2)))) + 5.4;
                Pshow = Pdb;  fShow = f;  % first block, for the plot
            end
            [a, c, w, e] = frameCarrier(Pdb, noiseDB, f, S);
            act = [act; a]; car = [car; c]; bwF = [bwF; w]; exF = [exF; e]; %#ok<AGROW>
        end
        nFr = numel(act);

        % --- B) Bursts, channels, hop groups ---
        Bu = segmentBursts(act, car, bwF, exF, S);
        if height(Bu) == 0
            fclose(fid);
            warning('  -> %s icinde burst bulunamadi, atlaniyor.', nm);
            continue;
        end

        Bu.Slice     = repmat(string(sliceTag), height(Bu), 1);
        Bu.tStart_s  = tOffset + S.tStart + (Bu.sFrame - 1) * S.win / Fs;
        Bu.dur_s     = Bu.nFrames * S.win / Fs;
        Bu.atStart   = Bu.sFrame == 1;
        Bu.atEnd     = Bu.eFrame == nFr;
        Bu.fcAbs_GHz = (fc + Bu.fcRel_Hz) / 1e9;

        minFrHop = ceil(S.minBurst_us*1e-6*Fs / S.win) + 2;
        isLong = Bu.nFrames >= minFrHop & Bu.exc_dB >= S.hopMinExc_dB;
        Bu = Bu(isLong, :);
        nB = height(Bu);

        if nB == 0
            fclose(fid);
            warning('  -> %s icinde uzun burst yok, atlaniyor.', nm);
            continue;
        end

        [ids, centers, counts] = assignChannels(Bu.fcRel_Hz, S.chanTol);
        Bu.chan = ids;
        gs = find([true; diff(ids) ~= 0]);
        ge = [gs(2:end) - 1; nB];
        Gp = table();
        for k = 1:numel(gs)
            t0 = Bu.tStart_s(gs(k));
            t1 = Bu.tStart_s(ge(k)) + Bu.dur_s(ge(k));
            Gp = [Gp; table(k, ids(gs(k)), (fc + centers(ids(gs(k))))/1e9, t0, t1, t1 - t0, ge(k) - gs(k) + 1, ...
                Bu.atStart(gs(k)), Bu.atEnd(ge(k)), ...
                'VariableNames', {'Group','Chan','Chan_GHz','t0_s','t1_s','Dwell_s','nBursts','atStart','atEnd'})]; %#ok<AGROW>
        end

        H = hopSummary(Bu, Gp, centers, counts, act);
        H_tbl = struct2table(H, 'AsArray', true);
        H_tbl.Slice = string(sliceTag);
        H_tbl = movevars(H_tbl, 'Slice', 'Before', 1);
        All_H = [All_H; H_tbl]; %#ok<AGROW>

        % --- C) Per-burst feature extraction ---
        minFr = ceil(S.minBurst_us*1e-6*Fs / S.win) + 2;
        sel = find(Bu.nFrames >= minFr & ~Bu.atStart & ~Bu.atEnd);
        if numel(sel) > S.maxBurstFeat
            sel = sel(unique(round(linspace(1, numel(sel), S.maxBurstFeat))));
        end

        Ft = table();
        for q = sel(:).'
            i1  = n0 + Bu.sFrame(q) * S.win;
            len = (Bu.nFrames(q) - 2) * S.win;
            x   = readCF32(fid, i1, len);
            xb  = x .* exp(-1j * 2*pi * Bu.fcRel_Hz(q) / Fs * (0:len-1).');
            y   = upfirdn(xb, hLP, 1, S.decim);
            try
                Fb = extract_rf_features(y, Fs2, fc + Bu.fcRel_Hz(q), 'trim', trimN, 'winSTFT', 256, 'usableBW', 18e6);
                Fb = rmfield(Fb, S.dropFields);
                Ft = [Ft; [table(string(sliceTag), q, Bu.tStart_s(q), Bu.dur_s(q), Bu.fcAbs_GHz(q), ...
                    'VariableNames', {'Slice','Burst','tStart_s','dur_s','fcAbs_GHz'}), struct2table(Fb)]]; %#ok<AGROW>
            catch ME
                warning('  -> %s Burst %d ozellik hatasi: %s', sliceTag, q, ME.message);
            end
        end

        % --- D) 4-panel overview for every 1 s slice ---
        fSlice = figure('Color', 'w', 'Position', [50 50 1400 900], 'Visible', 'off');
        try, theme(fSlice, 'light'); catch, end % R2025a+: avoid dark-theme (light-gray) text on the white figure
        tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

        nexttile;
        nShow = size(Pshow, 2);
        imagesc(tOffset + S.tStart + (0:nShow-1)*S.win/Fs*1e3, (fc + fShow)/1e9, Pshow);
        axis xy; colormap('jet');
        pk = max(Pshow(:)); caxis([pk - 60, pk]); colorbar;
        xlabel('Zaman (ms)'); ylabel('Frekans (GHz)');
        title(sprintf('Spektrogram (%s - ilk blok)', sliceTag));

        nexttile;
        tF = tOffset + S.tStart + ((0:nFr-1)' + 0.5) * S.win / Fs;
        scatter(tF(act)*1e3, (fc + car(act))/1e9, 6, exF(act), 'filled'); colorbar; grid on;
        xlabel('Zaman (ms)'); ylabel('Tasiyici (GHz)');
        title(sprintf('Kare basina tasiyici izi (%s)', sliceTag));

        nexttile;
        validDwell = Gp.Dwell_s(~Gp.atStart & ~Gp.atEnd)*1e3;
        if ~isempty(validDwell)
            histogram(validDwell, 20);
        end
        grid on;
        xlabel('Dwell (ms)'); ylabel('Hop sayisi');
        title('Kanalda kalis suresi (kenara degmeyen gruplar)');

        nexttile;
        bar((fc + centers)/1e9, counts); grid on;
        xlabel('Kanal (GHz)'); ylabel('Burst sayisi');
        title(sprintf('Kanal kullanimi (%d kanal)', numel(centers)));

        try, sgtitle(sprintf('%s | %s | %s', meta.Drone, char(rows.relFolder(1)), nm), 'Interpreter', 'none'); catch, end
        saveas(fSlice, fullfile(S.outDir, sprintf('overview_%s.png', sliceTag)));
        close(fSlice);

        fclose(fid);

        All_Ft = [All_Ft; Ft]; %#ok<AGROW>
        All_Bu = [All_Bu; Bu]; %#ok<AGROW>
        fprintf('  -> Tamamlandi: %d uzun burst, %d ozellik, grafik: overview_%s.png\n', nB, height(Ft), sliceTag);

        if S.savePerSlice && ~isempty(Ft)
            subOut = fullfile(S.outDir, sliceTag);
            if ~exist(subOut, 'dir'), mkdir(subOut); end
            writetable(Ft, fullfile(subOut, 'burst_features.csv'));
            writetable(H_tbl, fullfile(subOut, 'hop_summary.csv'));
        end
    end

    %% ---- 3) Mean and median over the whole recording ----
    assert(~isempty(All_Ft), 'Hicbir dosyadan gecerli ozellik cikarilamadi!');

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

    % Final feature vector in a single row
    Final_Combined_Vector = [Overall_Hop_Mean, Overall_Feat_Mean(:, 3:end)];

    %% ---- 4) Print and save ----
    fprintf('\n-----------------------------------------------------------------\n');
    fprintf('TUM DILIMLER (%d Dilim, %d Burst) GENEL ORTALAMA OZELLIKLER:\n', ...
        height(All_H), height(All_Ft));
    fprintf('-----------------------------------------------------------------\n');
    disp(Overall_Feat_Mean);
    disp(Overall_Hop_Mean);

    writetable(All_Bu,                fullfile(S.outDir, 'all_slices_bursts.csv'));
    writetable(All_H,                 fullfile(S.outDir, 'all_slices_hop_summary.csv'));
    writetable(All_Ft,                fullfile(S.outDir, 'all_slices_burst_features.csv'));
    writetable(Overall_Feat_Mean,     fullfile(S.outDir, 'overall_burst_features_mean.csv'));
    writetable(Overall_Feat_Median,   fullfile(S.outDir, 'overall_burst_features_median.csv'));
    writetable(Overall_Hop_Mean,      fullfile(S.outDir, 'overall_hop_summary_mean.csv'));
    writetable(Final_Combined_Vector, fullfile(S.outDir, 'FINAL_COMBINED_FEATURE_VECTOR.csv'));

    fprintf('Tum ciktilar kaydedildi: %s\n', S.outDir);

    res.nSlicesUsed = height(All_H);
    res.nBursts = height(All_Bu);
    res.nFeat = height(All_Ft);
    if ismember('snr_dB', All_Ft.Properties.VariableNames)
        res.snr_med = median(All_Ft.snr_dB, 'omitnan');
    end
    res.status = 'ok';
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

function x = readCF32(fid, startSample, n)
    fseek(fid, 8 * startSample, 'bof');
    raw = fread(fid, [2, n], 'float32=>double');
    x = complex(raw(1, :).', raw(2, :).');
end

function [act, car, bw, exc] = frameCarrier(Pdb, noiseDB, f, S)
    E = Pdb - noiseDB;
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