Barışın ve Rüzgar'ın erişimi olan, Yusifin açtığından farklı olan, bitirme projesi kodlarına erişeceğimiz repository. 
	Status
Hardware test bench	One ADALM-Pluto SDR (transmit and receive channels), antennas, host computer, USB connection	Available
Data and preprocessing	Analysis of recordings, conversion to SDR-compatible signals, SNR control	(planned)
SDR control and recording	Transmit/receive control, automatic recording and labeling	(planned)
Detection and classification	Energy/CFAR detection, feature extraction, SVM/Random Forest classifier	(planned)
Evaluation and demo	Metrics across SNR levels, end-to-end demonstration	(planned)
Performance targets

Measured under controlled laboratory conditions:

Metric	Target
Drone detection probability at SNR ≥ 0 dB	≥ 95 %
False alarm probability	≤ 1 %
Macro-F1 for basic classification at SNR ≥ 0 dB	≥ 85 %
Reported SNR range	−10 dB to +10 dB
Reproducibility	Same input recording and SDR settings give repeatable results
How SNR is controlled

To reach a target SNR, noise is added to the clean recording. The noise power is derived as follows.

Step 1: signal power. For a complex baseband signal $x[n]$ with $N$ samples, the average power is

$$P_s = \frac{1}{N}\sum_{n=0}^{N-1} |x[n]|^2$$

Step 2: SNR in linear units. The target SNR is given in dB, so convert it:

$$\mathrm{SNR}_{\mathrm{lin}} = 10^{\mathrm{SNR}_{\mathrm{dB}}/10}$$

Step 3: required noise power. By definition $\mathrm{SNR}_{\mathrm{lin}} = P_s / \sigma^2$, so

$$\sigma^2 = \frac{P_s}{\mathrm{SNR}_{\mathrm{lin}}}$$

Step 4: generate the noise. Complex white Gaussian noise $w[n] \sim \mathcal{CN}(0,\sigma^2)$ splits its power equally between the I and Q parts, so each part has variance $\sigma^2/2$:

$$w[n] = \sqrt{\frac{\sigma^2}{2}}\,\bigl(a[n] + j\,b[n]\bigr), \qquad a[n],\, b[n] \sim \mathcal{N}(0,1)$$

Step 5: noisy signal.

$$y[n] = x[n] + w[n]$$

This noise is added in software before transmission. Over-the-air or cable transmission then adds the real hardware behavior of the SDR chain.

Hardware and limits
SDR: ADALM-Pluto (one device; transmit and receive channels are used together).
Frequency range: stock Pluto covers about 325 MHz to 3.8 GHz, so only recordings in this range can be replayed. Recordings near 5.8 GHz are outside the stock range.
Bandwidth and sample rate: the Pluto's bandwidth is smaller than that of many recording instruments, so band selection and resampling are required. Check the current Analog Devices documentation for exact limits.
ATEL-provided items: SDR access, RF cables, attenuators and test environment, as agreed for the project (confirm what you actually received).
Safety and regulations
Transmit only at frequencies and power levels that you are permitted to use.
Prefer a wired connection with an attenuator between the transmit and receive ports. If you must transmit over the air, use low power, and check the national rules first.
Never connect the transmit port directly to the receive port without attenuation, since this can damage the receiver.
Use this project only for laboratory testing of detection algorithms.
