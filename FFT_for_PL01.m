close all 
clear all
clc

% fast forier transform 
% Ben Lowin
% 11/26/2024

%%
load PL01_NCP_02.mat

%% Define base varables    

T = 10;                              % Sampling period  (sec) 
Fs = 1/T;                            % Sampling frequency (Hz)
L = length(pigi_dat.datetime);       % Length of signal
t = (0:L-1)*T;                       % Time vector


%% need to fill missing data, using linerar interpolation

varable = 'NCP';

figure
plot(pigi_dat.datetime, pigi_dat.NCP)

inter_O2 = fillmissing(pigi_dat.NCP', 'linear'); % Linearly interpolate NaN values

%% Use the FFT

Y = fft (inter_O2);

powerspectrum_Y = abs(Y/L).^2;%

f = Fs*(0:(L/2))/L;


%% make a 95% threshold noise

shift_rand = mean(inter_O2);

% Noise model for significance threshold
noise = randn(size(inter_O2)) * shift_rand; % Generate random noise
fft_noise = fft(noise); % FFT of noise
power_spectrum_noise = abs(fft_noise/L).^2; % Power spectrum of noise

% Compute 95% significance threshold
threshold = prctile(power_spectrum_noise, 95);

%% Plot with threshold

%note that the current shows way to far
% one year  is 3.17*10^-8 Hz
% One Day   is 1.157*10^-5 Hz
% six hours is 4.6296e-05 Hz
% One Hour  is 2.77*10^-4 Hz

figure;
plot(f(1:floor(L/2)), powerspectrum_Y(1:floor(L/2)), 'b', 'LineWidth', 1.5); % Power spectrum
hold on;
yline(threshold, 'r--', 'LineWidth', 1.5); % Significance threshold
xlabel('Frequency (Hz)');
ylabel('Power');
title(['PL-01 - FFT ' varable]);
legend('Power Spectrum', '95% Threshold');
ylim([0 1000])
xlim([0 1*10^-4 ])

text(4.45815*10^-5, 85, '6.23 hours')
text(2.22907*10^-5, 380, '12.46 hours')
text(1.11454*10^-5, 522, '24.92 hours')


% Note that the M2 tide has a period of approximately 12.42 hours
peak_1 = 1/(4.45815*10^-5)/3600;
peak_2 = 1/(2.22907*10^-5)/3600;
peak_3 = 1/(1.11454*10^-5)/3600;




