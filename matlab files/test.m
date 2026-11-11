folder = 'C:\Users\User\Desktop\ATEL(1)\DJI_Inspire1Pro';
files  = dir(fullfile(folder,'*.mat'));     % only .mat files
K      = numel(files);   fprintf('%d files found\n',K);
fc     = nan(K,1);

for k = 1:K
    try
        d  = load(fullfile(files(k).folder,files(k).name)).data;
        fs = 1/d.XInc;                              % from the file, not hard-coded
        y  = double(d.Data(:))*d.YInc + d.YOrg;     % volts

        % burst start: first 1-us window with power above 10% of the maximum
        w  = round(1e-6*fs);
        pw = movmean(y.^2, w);
        n0 = find(pw > 0.1*max(pw), 1);

        % 5 us after the start, up to about 110 us, power-of-2 length
        seg = y(n0+round(5e-6*fs):end);
        N   = 2^floor(log2(min(numel(seg), round(110e-6*fs))));
        seg = seg(1:N) - mean(seg(1:N));

        X = abs(fft(seg.*hann(N))).^2;
        f = (0:N/2)'*fs/N;   X = X(1:N/2+1);
        m = f>2.3e9 & f<2.55e9;                     % widen if needed
        fc(k) = sum(f(m).*X(m)) / sum(X(m));

        % downconvert to 0 Hz and decimate to 20 MS/s
        z = y .* exp(-1j*2*pi*fc(k)*(0:numel(y)-1)'/fs);
        for p = factor(round(fs/20e6))
            z = decimate(real(z),p) + 1j*decimate(imag(z),p);
        end
        if k==1 || ~exist('IQ','var'), IQ = zeros(numel(z),K,'single'); end
        IQ(:,k) = single(z);
        fprintf('%d/%d  fc = %.4f GHz\n',k,K,fc(k)/1e9);
    catch ME
        warning('skipped %s: %s', files(k).name, ME.message);
    end
end

figure; plot(fc/1e6,'.-'); xlabel('reading #'); ylabel('carrier (MHz)'); grid on

%% spectrogram

fsb = 20e6;                                  % baseband rate of IQ
fg  = (2.40e9 : 78.125e3 : 2.48e9)';          % common RF grid, 78 kHz steps
M   = zeros(numel(fg), K);

for k = 1:K
    if isnan(fc(k)), continue, end
    [P,f] = pwelch(IQ(2600:4900,k), hann(256), 128, 256, fsb, 'centered');
    M(:,k) = interp1(fc(k)+f, P, fg, 'linear', 0);   % place on absolute RF axis
end

imagesc(1:K, fg/1e6, 10*log10(M + eps));
axis xy; colormap turbo; colorbar
xlabel('reading #'); ylabel('RF frequency (MHz)'); title('Hop map');
%% 
% This shows which channels the drone uses and whether use is uniform, which is a key feature for classification (WP5).
edges = 2400:2:2480;
histogram(fc/1e6, edges); xlabel('carrier (MHz)'); ylabel('# readings');