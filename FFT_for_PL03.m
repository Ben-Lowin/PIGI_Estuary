close all 
clear all
clc

% fast forier transform 
% Ben Lowin
% 11/26/2024

%%
load PL03_NCP_01.mat

% there are missing time frames in this data set
% this needs to be corrected before the FFT can be used

%% Correcting for NaN, NaT and time skips
var = 'NCP';

%start to end in one min intervals
continuios_time = pigi_dat.datetime(350):seconds(10):pigi_dat.datetime(end-14);

% O2 values have NaN, need to remove
% there are also 14 not a time values, these need to be removed as well
input  = pigi_dat.NCP(350:end-14);
%  input(52705:52731) = NaN; % for salinity

inter_O2 = fillmissing(input, 'linear'); % Linearly interpolate NaN values

time_p = pigi_dat.datetime(350:end-14); %with out the last 14 NaT values


%interpolate the O2 data to this, will also fill NaNs?
continuios_O2 = interp1(time_p, inter_O2, continuios_time);


figure
plot (continuios_time, continuios_O2)

%% Define base varables    

T = 10;                              % Sampling period  (sec) 
Fs = 1/T;                            % Sampling frequency (Hz)
L = length(continuios_time);         % Length of signal
t = (0:L-1)*T;                       % Time vector



%% Use the FFT

Y = fft (continuios_O2);
powerspectrum_Y = abs(Y/L).^2;%

f = Fs*(0:(L/2))/L; % frequency


%% make a 95% threshold noise

shift_rand = mean(continuios_O2);

% Noise model for significance threshold
noise = randn(size(continuios_O2)) * shift_rand; % Generate random noise
fft_noise = fft(noise); % FFT of noise
power_spectrum_noise = abs(fft_noise/L).^2; % Power spectrum of noise

% Compute 95% significance threshold
threshold = prctile(power_spectrum_noise, 95);

%% Plot with threshold

%note that the current shows way to far
% one year in Hz is 3.17*10^-8 Hz
% One Day in Hz is 1.157*10^-5 Hz
% six hours is 4.6296e-05 Hz
% One Hour in Hz is 2.77*10^-4 Hz

figure;
plot(f(1:floor(L/2)), powerspectrum_Y(1:floor(L/2)), 'b', 'LineWidth', 1.5); % Power spectrum
hold on;
yline(threshold, 'r--', 'LineWidth', 1.5); % Significance threshold
xlabel('Frequency (Hz)');
ylabel('Power');
title(['PL-03 - FFT ' var]);
legend('Power Spectrum', '95% Threshold');
% ylim([0 2])
xlim([0 20*10^-5 ])

text(1.13235*10^-5, 8280, '24.53 hours')
text(2.26471*10^-5, 2750, '12.26 hours')
text(2.2647*10^-5, .6, '12.26 hours')



% Note that the M2 tide has a period of approximately 12.42 hours
peak_1 = 1/(1.13235*10^-5)/3600;
peak_2 = 1/(2.26471*10^-5)/3600;
peak_3 = 1/(2.2647*10^-5)/3600;
peak_4 = 1/(5.66177*10^-6)/3600;
peak_5 = 1/(9.43628*10^-6)/3600;




