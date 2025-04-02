close all 
clear all
clc

% Fast Fourier Transform with Spectrogram and Cross-Spectra
% Ben Lowin
% Updated: 03/13/2025

%% Load Data
load PL01_NCP_03.mat

%% Define Base Variables    
T = 10;                              % Sampling period  (sec) 
Fs = 1/T;                            % Sampling frequency (Hz)
L = length(pigi_dat.datetime);       % Length of signal
t = (0:L-1)*T;                       % Time vector

%% Fill Missing Data Using Linear Interpolation
var_name = 'salintiy';
figure;
plot(pigi_dat.datetime, pigi_dat.ts_cor.sal);
title('Original NCP Data');
xlabel('Time');
ylabel('Salinity');

inter_O2 = fillmissing(pigi_dat.ts_cor.sal', 'linear'); % Linearly interpolate NaN values

%% Apply FFT with Windowing (Hanning)
window = hann(L)'; % Hanning Window
Y = fft(inter_O2 .* window); % Apply windowing before FFT

powerspectrum_Y = abs(Y/L).^2; % Power Spectrum
f = Fs*(0:(L/2))/L; % Frequency Vector

%% Compute Noise Threshold (95%)
shift_rand = mean(inter_O2);
noise = randn(size(inter_O2)) * shift_rand; % Random Noise
fft_noise = fft(noise .* window); % Apply windowing to noise
power_spectrum_noise = abs(fft_noise/L).^2; % Noise Power Spectrum
threshold = prctile(power_spectrum_noise, 95); % 95% threshold

%% Plot FFT with 95% Threshold
figure;
plot(f(1:floor(L/2)), powerspectrum_Y(1:floor(L/2)), 'b', 'LineWidth', 1.5);
hold on;
yline(threshold, 'r--', 'LineWidth', 1.5); % Significance threshold
xlabel('Frequency (Hz)');
ylabel('Power');
title(['PL-01 - FFT ', var_name, ' with Hanning Window']);
legend('Power Spectrum', '95% Threshold');
ylim([0 1000])
xlim([0 1*10^-4 ])

text(4.45815*10^-5, 85, '6.23 hours')
text(2.22907*10^-5, 380, '12.46 hours')
text(1.11454*10^-5, 522, '24.92 hours')

%% Spectrogram Analysis (Time-Frequency)
window_size = round(L / 3); % Define Window Size (~10 segments)
overlap = round(window_size * 0.5); % 50% Overlap
nfft = 2^nextpow2(window_size); % Optimized FFT Size

figure;
spectrogram(inter_O2, hann(window_size), overlap, nfft, Fs, 'yaxis');
title('Spectrogram of NCP');
ylabel('Frequency (Hz)');
xlabel('Time');
colorbar;

%% Cross-Spectral Analysis (Covariance at Certain Frequencies)
if isfield(pigi_dat, 'ts_cor') % Replace with a real variable
    other_variable = fillmissing(pigi_dat.ts_cor.solar_rad', 'linear'); % Interpolate missing data
    
    % Compute cross-power spectral density
    [Pxy, F] = cpsd(inter_O2, other_variable, hann(window_size), overlap, nfft, Fs);
    
    % Plot Cross-Spectrum
    figure;
    plot(F, abs(Pxy), 'LineWidth', 1.5);
    title('Cross-Spectral Density: NCP vs. Another Variable');
    xlabel('Frequency (Hz)');
   
    ylabel('Power');
    grid on;
end
