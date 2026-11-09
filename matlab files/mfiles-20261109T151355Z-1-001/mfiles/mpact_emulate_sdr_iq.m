function info = mpact_emulate_sdr_iq(matFile, outIqFile, fc, Fs_sdr, bw, opts)
%MPACT_EMULATE_SDR_IQ  MPACT DroneRC RF dataset real waveform -> SDR-like complex IQ emulation.
%
%   info = mpact_emulate_sdr_iq(matFile, outIqFile, fc, Fs_sdr, bw, opts)
%
% INPUTS
%   matFile   : MPACT .mat file path, containing 'data' struct (data.Data, data.YInc, data.NumPoints, ...)
%   outIqFile : output file path for SDR IQ (sc16 int16 interleaved: I0,Q0,I1,Q1,...)
%   fc        : desired downconversion center (Hz). For this dataset typically ~2.4e9 (or band center).
%   Fs_sdr    : target SDR sampling rate (Hz), e.g., 20e6
%   bw        : target channel bandwidth (Hz), e.g., 20e6 (LPF cutoff will be based on this)
%   opts      : options struct (see defaults below)
%
% OUTPUT
%   info struct with metadata and achieved lengths.
%
% IMPORTANT NOTES
% - Dataset is real (single channel). We create an analytic signal using Hilbert transform.
% - If your MAT file has a field data.XInc or data.SamplingRate, set opts.Fs_raw from it. Default is 20e9.
% - For huge files, enable cropping (opts.crop=true) to process manageable segments.
%
% Requires: Signal Processing Toolbox for hilbert(), resample(), fir1().

arguments
    matFile (1,:) char
    outIqFile (1,:) char
    fc (1,1) double
    Fs_sdr (1,1) double {mustBePositive}
    bw (1,1) double {mustBePositive}

    opts.Fs_raw (1,1) double {mustBePositive} = 20e9      % dataset doc: 20 GSa/s
    opts.useVolts (1,1) logical = false                   % true: x <- x * data.YInc
    opts.removeMean (1,1) logical = true                  % DC removal
    opts.crop (1,1) logical = true
    opts.cropStart (1,1) double = 1                       % sample index (1-based)
    opts.cropLength (1,1) double = 2e6                    % number of samples to process

    % Optional pre-bandpass before Hilbert (recommended if you know band)
    opts.preBandpass (1,1) logical = false
    opts.preBP_bw (1,1) double {mustBePositive} = 200e6    % bandpass BW around fc at Fs_raw

    % LPF after DDC (anti-alias)
    opts.lpfTaps (1,1) double {mustBeInteger, mustBePositive} = 4096
    opts.lpfRollOff (1,1) double {mustBePositive} = 0.8    % cutoff = rollOff*(bw/2)

    % Impairments
    opts.addDC (1,1) logical = false
    opts.dcI (1,1) double = 0.001
    opts.dcQ (1,1) double = -0.0007

    opts.addIQImbalance (1,1) logical = false
    opts.imageLeakAlpha (1,1) double = 0.01               % y = y + alpha*conj(y)

    opts.addPhaseNoise (1,1) logical = false
    opts.phaseNoiseStd (1,1) double = 5e-5                % rad/sample random-walk increment

    % Quantization / export
    opts.exportFormat (1,:) char = 'sc16'                 % 'sc16' or 'cf32'
    opts.headroom (1,1) double {mustBePositive} = 0.8     % scale to this max magnitude before quantize
end

% ----------------- Load MAT -----------------
S = load(matFile, 'data');
assert(isfield(S,'data'), "MAT dosyasında 'data' struct bulunamadı.");
data = S.data;

assert(isfield(data,'Data'), "data.Data bulunamadı.");
x = double(data.Data(:));

% If available, prefer real Fs from MAT metadata
if isfield(data,'XInc') && ~isempty(data.XInc) && data.XInc > 0
    opts.Fs_raw = 1.0 / double(data.XInc);
end
if isfield(data,'SamplingRate') && ~isempty(data.SamplingRate) && data.SamplingRate > 0
    opts.Fs_raw = double(data.SamplingRate);
end

Fs_raw = opts.Fs_raw;

% ----------------- Preprocess -----------------
if opts.removeMean
    x = x - mean(x);
end

if opts.useVolts
    assert(isfield(data,'YInc'), "opts.useVolts=true ama data.YInc yok.");
    x = x * double(data.YInc);
end

% Crop
N0 = numel(x);
if opts.crop
    s0 = max(1, round(opts.cropStart));
    L  = round(opts.cropLength);
    s1 = min(N0, s0 + L - 1);
    x = x(s0:s1);
else
    s0 = 1;
    s1 = N0;
end

% ----------------- Optional pre-bandpass around fc -----------------
% This improves analytic signal quality if the raw capture includes a lot of out-of-band content.
if opts.preBandpass
    bp_half = opts.preBP_bw/2;
    f1 = max(1, fc - bp_half);
    f2 = min((Fs_raw/2)*0.999, fc + bp_half);
    % Design real bandpass FIR at Fs_raw
    Wn = [f1 f2] / (Fs_raw/2);
    bpTaps = min(8192, max(1024, 2*ceil(opts.lpfTaps/2))); % reasonable order
    hbp = fir1(bpTaps-1, Wn, 'bandpass', kaiser(bpTaps,7));
    x = filtfilt(hbp, 1, x);  % zero-phase for analysis
end

% ----------------- Real -> Analytic (complex) -----------------
xa = hilbert(x);  % analytic signal

% ----------------- DDC -----------------
n = (0:numel(xa)-1).';
lo = exp(-1j*2*pi*fc/Fs_raw * n);
xbb = xa .* lo;

% ----------------- LPF (anti-alias / channel select) -----------------
cutoff = opts.lpfRollOff * (bw/2);
cutoff = min(cutoff, 0.49*(Fs_raw/2));
wn = cutoff/(Fs_raw/2);

h = fir1(opts.lpfTaps-1, wn, 'low', kaiser(opts.lpfTaps, 7));
xbb_f = filter(h, 1, xbb);

% ----------------- Resample -----------------
[p,q] = rat(Fs_sdr/Fs_raw, 1e-12);
y = resample(xbb_f, p, q);

% ----------------- Impairments (optional) -----------------
if opts.addDC
    y = y + complex(opts.dcI, opts.dcQ);
end

if opts.addIQImbalance
    y = y + opts.imageLeakAlpha * conj(y);
end

if opts.addPhaseNoise
    dphi = opts.phaseNoiseStd * randn(size(y));
    phi = cumsum(dphi);
    y = y .* exp(1j*phi);
end

% ----------------- Export -----------------
switch lower(opts.exportFormat)
    case 'sc16'
        % Normalize to headroom
        m = max(abs(y));
        if m > 0
            y = y * (opts.headroom / m);
        end
        I = int16(max(min(round(real(y)*32767), 32767), -32768));
        Q = int16(max(min(round(imag(y)*32767), 32767), -32768));

        iq = zeros(2*numel(I),1,'int16');
        iq(1:2:end) = I;
        iq(2:2:end) = Q;

        fid = fopen(outIqFile,'wb');
        assert(fid>0, "Çıktı dosyası açılamadı: %s", outIqFile);
        fwrite(fid, iq, 'int16');
        fclose(fid);

    case 'cf32'
        % interleaved float32 (I,Q)
        tmp = [single(real(y)) single(imag(y))].';
        fid = fopen(outIqFile,'wb');
        assert(fid>0, "Çıktı dosyası açılamadı: %s", outIqFile);
        fwrite(fid, tmp(:), 'float32');
        fclose(fid);

    otherwise
        error("Bilinmeyen exportFormat: %s (sc16/cf32)", opts.exportFormat);
end

% ----------------- Return info -----------------
info = struct();
info.matFile = matFile;
info.outIqFile = outIqFile;
info.Fs_raw = Fs_raw;
info.fc = fc;
info.Fs_sdr = Fs_sdr;
info.bw = bw;
info.inputSamples = numel(x);
info.inputSampleRange = [s0 s1];
info.outputSamples = numel(y);
info.resample_pq = [p q];
info.exportFormat = opts.exportFormat;

fprintf("Done.\n");
fprintf("  Fs_raw = %.6f GHz, fc = %.6f GHz\n", Fs_raw/1e9, fc/1e9);
fprintf("  Input samples (used) = %d (range %d..%d)\n", info.inputSamples, s0, s1);
fprintf("  Fs_sdr = %.3f MHz, Output samples = %d\n", Fs_sdr/1e6, info.outputSamples);
fprintf("  Export = %s -> %s\n", opts.exportFormat, outIqFile);

end
