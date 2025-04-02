close all
clear all
clc

% Figures for the Three Methods Paper
% Ben Lowin
% Febuary 11th 2025

%%
% Define common axis formatting
fontSize = 12;
fontName = 'Times New Roman';
lineWidth = 1.2;


%% Figure 1: Map of PL with Georgia popout
%struggling with this map - moving on 
% https://www.ncei.noaa.gov/products/etopo-global-relief-model
% Extract a grid from the NCEI Coastal Relief Models (CRM) Mosaic,
% a bathymetry/topography mosaic that includes NCEI-stewarded DEMs. 
% The depth values are in meters, stored as 32-bit floating point values. 
% The cell size of the CRM is 1 arc-second. More information about Coastal Relief Models at NCEI

[depth_eto,spatial_ref_info] = readgeoraster('Wassaw_sound_map.tif');%reads in data
[rows,cols] = size(depth_eto);%pulls the size
[rowGrid, colGrid] = ndgrid(1:rows, 1:cols); %replicates the grid
[YY_eto, XX_eto] = intrinsicToGeographic(spatial_ref_info, colGrid, rowGrid);%changes to lat lon
depth_eto = double(depth_eto);%formats to double

moutain = 1000;
depth_min = 0 ;
depth_max = -5 ;


ax1 = figure;
[C,h] = contourf(XX_eto, YY_eto, depth_eto,[moutain,depth_min]) ;
hold on
text(-81.0421, 31.9438, 'Preists Landing','FontSize',12,'FontName',...
    'Times New Roman','Rotation',45)
plot(-81.013708, 31.963728,  'ro', 'MarkerSize',10,'MarkerFaceColor','r')
text(-80.9679, 31.9313, sprintf('Wassaw \n Sound'),'FontSize',16,'FontName',...
    'Times New Roman')
text(-80.9679, 31.9313, sprintf('Wassaw \n Sound'),'FontSize',16,'FontName',...
    'Times New Roman')

hold off
colormap([0.7, 0.9, 0.7]);
xlabel('Longitude');
ylabel('Latitude');


fig = gcf;
fig.Position = [100, 100, 1000, 800];


% print(fig, '-dpng', '-r600', 'high_res_map.png'); % PNG (for self)
% print(fig, 'high_res_map', '-depsc', '-r600'); % EPS (for publication)

clear all
%% Figure 2 - Time Series for the June Data (Publication Quality)

% Define common axis formatting
fontSize = 12;
fontName = 'Times New Roman';
lineWidth = 1.2;

load PL01_NCP_03.mat

% Extract variables
oxygen  = pigi_dat.qc_cal.o2uM;
tide = pigi_dat.tide_depth;
NCP = pigi_dat.NCP ./ mean(tide, 'omitnan');
sal = pigi_dat.qc_cal.sal;
par = pigi_dat.ts_cor.solar_rad;
time = pigi_dat.datetime;

% Create figure with optimized layout
figure
tiledlayout(4,1, 'TileSpacing', 'compact', 'Padding', 'compact'); % Minimize spacing



% First subplot
ax1 = nexttile;
plot(time, oxygen, 'LineWidth', lineWidth)
hold on
ylabel('O_2 (mM)', 'FontSize', fontSize, 'FontName', fontName)
set(gca, 'XTickLabel', [], 'FontSize', fontSize, 'FontName', fontName)
text(time(end-15),200,'A','FontSize', fontSize, 'FontName', fontName)
hold off

% Second subplot
ax2 = nexttile;
hold on
plot(time, NCP, 'LineWidth', lineWidth)
ylabel('NCP (mmolO_2/m^2/d)', 'FontSize', fontSize, 'FontName', fontName)
set(gca, 'XTickLabel', [], 'FontSize', fontSize, 'FontName', fontName)
text(time(end-15),-25,'B','FontSize', fontSize, 'FontName', fontName)
hold off

% Third subplot
ax3 = nexttile;
yyaxis left
plot(time, sal, 'LineWidth', lineWidth)
ylabel('Sal (g/kg)', 'FontSize', fontSize, 'FontName', fontName)
text(time(end-15),29.5,'C','FontSize', fontSize, 'FontName', fontName)

yyaxis right
plot(time, tide, 'LineWidth', lineWidth)
ylabel('Water depth (m)', 'FontSize', fontSize, 'FontName', fontName)
set(gca, 'XTickLabel', [], 'FontSize', fontSize, 'FontName', fontName)

% Fourth subplot
ax4 = nexttile;
plot(time, par, 'LineWidth', lineWidth)
ylabel('Surface PAR (W/m^2)', 'FontSize', fontSize, 'FontName', fontName)
xlabel('Time', 'FontSize', fontSize, 'FontName', fontName) % Only bottom plot has x-label
text(time(end-15),1200,'D','FontSize', fontSize, 'FontName', fontName)

% Link x-axes
linkaxes([ax1, ax2, ax3, ax4], 'x');

% Adjust figure size for publication
fig = gcf;
fig.Position = [100, 100, 1000, 800];

% Ensure high-quality rendering
set([ax1, ax2, ax3, ax4], 'FontSize', fontSize, 'FontName', fontName, 'LineWidth', lineWidth, 'Box', 'on')

% print(fig, '-dpng', '-r600', 'June_time_series_02.png'); % PNG (for self)
% print(fig, 'June_time_series_02', '-depsc', '-r600'); % EPS (for publication)

% clear ax1 ax2 ax3 ax4 fig NCP oxygen par pigi_dat sal tide time
%% Figure 3 - Time Series for the June Data (Publication Quality)

load PL02_NCP_03.mat

% Extract variables
oxygen  = pigi_dat.qc_cal.o2uM;
tide = pigi_dat.tide_depth;
NCP = pigi_dat.NCP ./ mean(tide, 'omitnan');
sal = pigi_dat.qc_cal.sal;
par = pigi_dat.ts_cor.solar_rad;
time = pigi_dat.datetime;

% Create figure with optimized layout
figure
tiledlayout(4,1, 'TileSpacing', 'compact', 'Padding', 'compact'); % Minimize spacing



% First subplot
ax1 = nexttile;
plot(time, oxygen, 'LineWidth', lineWidth)
hold on
ylabel('O_2 (mM)', 'FontSize', fontSize, 'FontName', fontName)
set(gca, 'XTickLabel', [], 'FontSize', fontSize, 'FontName', fontName)
text(time(end-15),190,'A','FontSize', fontSize, 'FontName', fontName)
hold off

% Second subplot

xBox = [ time(53312), time(53312), time(end-15), time(end-15) ];
yBox = [-300 0 0 -300];

ax2 = nexttile;
fill(xBox, yBox, [0.8 0.8 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'none'); % RGB [0.5 0.5 0.5] = grey

hold on
plot(time, NCP, 'Color',[0 0.4470 0.7410], 'LineWidth', lineWidth)
ylabel('NCP (mmolO_2/m^2/d)', 'FontSize', fontSize, 'FontName', fontName)
set(gca, 'XTickLabel', [], 'FontSize', fontSize, 'FontName', fontName)
text(time(end-15),-25,'B','FontSize', fontSize, 'FontName', fontName)
hold off

% Third subplot
ax3 = nexttile;
yyaxis left
plot(time, sal, 'LineWidth', lineWidth)
ylabel('Sal (g/kg)', 'FontSize', fontSize, 'FontName', fontName)
text(time(end-15),30,'C','FontSize', fontSize, 'FontName', fontName)

yyaxis right
plot(time, tide, 'LineWidth', lineWidth)
ylabel('Water Dpeth (m)', 'FontSize', fontSize, 'FontName', fontName)
set(gca, 'XTickLabel', [], 'FontSize', fontSize, 'FontName', fontName)

% Fourth subplot
ax4 = nexttile;
plot(time, par, 'LineWidth', lineWidth)
ylabel('Surface PAR (W/m^2)', 'FontSize', fontSize, 'FontName', fontName)
xlabel('Time', 'FontSize', fontSize, 'FontName', fontName) % Only bottom plot has x-label
text(time(end-15),1000,'D','FontSize', fontSize, 'FontName', fontName)

% Link x-axes
linkaxes([ax1, ax2, ax3, ax4], 'x');

% Adjust figure size for publication
fig = gcf;
fig.Position = [100, 100, 1000, 800];

% Ensure high-quality rendering
set([ax1, ax2, ax3, ax4], 'FontSize', fontSize, 'FontName', fontName, 'LineWidth', lineWidth, 'Box', 'on')

% print(fig, '-dpng', '-r600', 'september_time_series_02.png'); % PNG (for self)
% print(fig, 'september_time_series_02', '-depsc', '-r600'); % EPS (for publication)

% clear ax1 ax2 ax3 ax4 fig NCP oxygen par pigi_dat sal tide time xBox yBox

%% Figure 4 - Time Series for the June Data (Publication Quality)

load PL03_NCP_03.mat

% Extract variables
oxygen  = pigi_dat.qc_cal.o2uM;
tide = pigi_dat.tide_depth;
NCP = pigi_dat.NCP ./ mean(tide, 'omitnan');
sal = pigi_dat.qc_cal.sal;
par = pigi_dat.ts_cor.solar_rad;
time = pigi_dat.datetime;

% Create figure with optimized layout
figure
tiledlayout(4,1, 'TileSpacing', 'compact', 'Padding', 'compact'); % Minimize spacing

% Define common axis formatting
fontSize = 12;
fontName = 'Times New Roman';
lineWidth = 1.2;

% First subplot
ax1 = nexttile;
plot(time, oxygen, 'LineWidth', lineWidth)
hold on
ylabel('O_2 (mM)', 'FontSize', fontSize, 'FontName', fontName)
set(gca, 'XTickLabel', [], 'FontSize', fontSize, 'FontName', fontName)
text(time(end-15),275,'A','FontSize', fontSize, 'FontName', fontName)
hold off

% Second subplot
ax2 = nexttile;
hold on
plot(time, NCP, 'LineWidth', lineWidth)
ylabel('NCP (mmolO_2/m^2/d)', 'FontSize', fontSize, 'FontName', fontName)
set(gca, 'XTickLabel', [], 'FontSize', fontSize, 'FontName', fontName)
text(time(end-15),-25,'B','FontSize', fontSize, 'FontName', fontName)
hold off

% Third subplot
ax3 = nexttile;
yyaxis left
plot(time, sal, 'LineWidth', lineWidth)
ylabel('Sal (g/kg)', 'FontSize', fontSize, 'FontName', fontName)
text(time(end-15),28,'C','FontSize', fontSize, 'FontName', fontName)

yyaxis right
plot(time, tide, 'LineWidth', lineWidth)
ylabel('Water Depth (m)', 'FontSize', fontSize, 'FontName', fontName)
set(gca, 'XTickLabel', [], 'FontSize', fontSize, 'FontName', fontName)

% Fourth subplot
ax4 = nexttile;
plot(time, par, 'LineWidth', lineWidth)
ylabel('Surface PAR (W/m^2)', 'FontSize', fontSize, 'FontName', fontName)
xlabel('Time', 'FontSize', fontSize, 'FontName', fontName) % Only bottom plot has x-label
text(time(end-15),600,'D','FontSize', fontSize, 'FontName', fontName)

% Link x-axes
linkaxes([ax1, ax2, ax3, ax4], 'x');

% Adjust figure size for publication
fig = gcf;
fig.Position = [100, 100, 1000, 800];

% Ensure high-quality rendering
set([ax1, ax2, ax3, ax4], 'FontSize', fontSize, 'FontName', fontName, 'LineWidth', lineWidth, 'Box', 'on')

% print(fig, '-dpng', '-r600', 'december_time_series_02.png'); % PNG (for self)
% print(fig, 'december_time_series_02', '-depsc', '-r600'); % EPS (for publication)

% clear ax1 ax2 ax3 ax4 fig NCP oxygen par pigi_dat sal tide time
%% Figure 4
%loads in all three data sets
load PL01_NCP_03.mat
PL01 = pigi_dat;

load PL02_NCP_03.mat
PL02 = pigi_dat;

load PL03_NCP_03.mat
PL03 = pigi_dat;

load PL03_three_plot.mat
load PL02_three_plot.mat
load PL01_three_plot.mat

clear pigi_dat

%read in sunrise/sunset times
j_s = readtable("June_time_sun.xlsx");
june_time = datenum(datetime(j_s.JuneTimes));

s_s = readtable("september_time_sun.xlsx");
sept_time = datenum(datetime(s_s.SeptemberTimes));

d_s = readtable("december_time_sun.xlsx");
dec_time = datenum(datetime(d_s.DecemberTimes));

clear j_s s_s d_s


%
figure
tiledlayout(3,1, 'TileSpacing', 'compact', 'Padding', 'compact');

% First subplot
ax1 = nexttile;
yyaxis left
hold on
% Add light-colored background boxes
for i = 1:2:length(june_time)-5
    rectangle('Position', [june_time(i), -300, june_time(i+1) - june_time(i), 300], ...
              'FaceColor', [0.9, 0.9, 0.9, 0.5], 'EdgeColor', 'none'); % Light gray with transparency
end
% Plot NCP data
plot(datenum(PL01.datetime), PL01.NCP./mean(PL01.tide_depth,'omitnan'), 'LineWidth', lineWidth)
plot(datenum(PL_01_three.odum_time),PL_01_three.odum_ncp, 'Color',[0, 0, 0],...
    'LineStyle','-','Marker','o','LineWidth',1.2,'MarkerFaceColor','k')

ylabel('NCP (mmolO_2/m^2/d)', 'FontSize', fontSize, 'FontName', fontName)
ylim([-300 0])

yyaxis right
% Plot tide depth
plot(datenum(PL01.datetime), PL01.tide_depth, 'LineWidth', lineWidth)
ylabel('Water Depth (m)', 'FontSize', fontSize, 'FontName', fontName)
set(gca, 'FontSize', fontSize, 'FontName', fontName)
ylim([0 8])

text(june_time(2)-.1,9,'A', 'FontSize', fontSize, 'FontName', fontName)

% Format x-axis as datetime
datetick('x', 'mmm-dd', 'keeplimits')
xlim([june_time(1) june_time(15)])
hold off

% Second subplot
ax2 = nexttile;
yyaxis left
hold on
% Add light-colored background boxes
for i = 1:2:length(sept_time)-3
    rectangle('Position', [sept_time(i), -300, sept_time(i+1) - sept_time(i), 300], ...
              'FaceColor', [0.9, 0.9, 0.9, 0.5], 'EdgeColor', 'none'); % Light gray with transparency
end

plot(datenum(PL02.datetime), PL02.NCP./mean(PL02.tide_depth,'omitnan'), 'LineWidth', lineWidth)
plot(datenum(PL_02_three.odum_time),PL_02_three.odum_ncp, 'Color',[0, 0, 0],...
    'LineStyle','-','Marker','o','LineWidth',1.2,'MarkerFaceColor','k')
ylabel('NCP (mmolO_2/m^2/d)', 'FontSize', fontSize, 'FontName', fontName)
ylim([-300 0])

yyaxis right
plot(datenum(PL02.datetime), PL02.tide_depth, 'LineWidth', lineWidth)
ylabel('Water Depth (m)', 'FontSize', fontSize, 'FontName', fontName)
set(gca, 'FontSize', fontSize, 'FontName', fontName)
ylim([0 8])

% Convert x-axis ticks to date format
datetick('x', 'mmm-dd', 'keepticks', 'keeplimits')
xlim([sept_time(1) sept_time(25)])
text(sept_time(2),9,'B', 'FontSize', fontSize, 'FontName', fontName)
hold off


% third subplot
ax3 = nexttile;
yyaxis left
hold on

% Add light-colored background boxes
for i = 1:2:length(dec_time)-5
    rectangle('Position', [dec_time(i), -300, dec_time(i+1) - dec_time(i), 300], ...
              'FaceColor', [0.9, 0.9, 0.9, 0.5], 'EdgeColor', 'none'); % Light gray with transparency
end

plot(datenum(PL03.datetime), PL03.NCP./mean(PL03.tide_depth,'omitnan'), 'LineWidth', lineWidth)
plot(datenum(PL_03_three.odum_time),PL_03_three.odum_ncp, 'Color',[0, 0, 0],...
    'LineStyle','-','Marker','o','LineWidth',1.2,'MarkerFaceColor','k')
ylabel('NCP (mmolO_2/m^2/d)', 'FontSize', fontSize, 'FontName', fontName)
ylim([-300 0])

yyaxis right
plot(datenum(PL03.datetime), PL03.tide_depth, 'LineWidth', lineWidth)
ylabel('Water Depth (m)', 'FontSize', fontSize, 'FontName', fontName)
set(gca, 'FontSize', fontSize, 'FontName', fontName)
ylim([0 8])

% Convert x-axis ticks to date format
datetick('x', 'mmm-dd', 'keepticks', 'keeplimits')
xlim([dec_time(1) dec_time(15)])
text(dec_time(2)-.35,9,'C', 'FontSize', fontSize, 'FontName', fontName)
hold off

% Adjust figure size for publication
fig = gcf;
fig.Position = [100, 100, 1000, 800];

% Ensure high-quality rendering
set([ax1, ax2, ax3], 'FontSize', fontSize, 'FontName', fontName, 'LineWidth', lineWidth, 'Box', 'on')

% print(fig, '-dpng', '-r600', 'NCP_PIGI_02.png'); % PNG (for self)
% print(fig, 'NCP_PIGI_02', '-depsc', '-r600'); % EPS (for publication)

% clear ax1 ax2 ax3 fig dec_time sept_time june_time i

%% Figure 5

load PL03_three_plot.mat
load PL02_three_plot.mat
load PL01_three_plot.mat

figure;
tiledlayout(1,3, 'TileSpacing', 'tight', 'Padding', 'compact');

%tile one
ax1 = nexttile;
errorbar(PL_01_three.btl_time, PL_01_three.btl_ncp, PL_01_three.btl_std)
hold on
errorbar( PL_01_three.pigi_time,PL_01_three.pigi_ncp,PL_01_three.pigi_std)
errorbar(PL_01_three.odum_time,PL_01_three.odum_ncp,PL_01_three.odum_std,'k')
text(PL_01_three.btl_time(9),-110,'A','FontSize', fontSize, 'FontName', fontName)

% legend('Bottle NCP', 'PIGI NCP','Odum NCP','FontSize', fontSize, 'FontName', fontName)
% xlabel('Date','FontSize', fontSize, 'FontName', fontName)
ylabel('NCP (mmolO_2/m^2/day)','FontSize', fontSize, 'FontName', fontName)
set(gca, 'FontSize', fontSize, 'FontName', fontName)
ylim([-125 75])
hold off


%tile one
ax2 = nexttile;
errorbar(PL_02_three.btl_time, PL_02_three.btl_ncp, PL_02_three.btl_std)
hold on
errorbar( PL_02_three.pigi_time,PL_02_three.pigi_ncp,PL_02_three.pigi_std)
errorbar(PL_02_three.odum_time,PL_02_three.odum_ncp,PL_02_three.odum_std,'k')
text(PL_02_three.btl_time(17),-110,'B','FontSize', fontSize, 'FontName', fontName)

% legend('Bottle NCP', 'PIGI NCP','Odum NCP','FontSize', fontSize, 'FontName', fontName)
% xlabel('Date','FontSize', fontSize, 'FontName', fontName)
% ylabel('NCP (mM/m^2/day)','FontSize', fontSize, 'FontName', fontName)
set(gca, 'YTickLabel', [], 'FontSize', fontSize, 'FontName', fontName)
ylim([-125 75])
hold off

%tile one
ax3 = nexttile;
errorbar(PL_03_three.btl_time, PL_03_three.btl_ncp, PL_03_three.btl_std)
hold on
errorbar( PL_03_three.pigi_time,PL_03_three.pigi_ncp,PL_03_three.pigi_std)
errorbar(PL_03_three.odum_time,PL_03_three.odum_ncp,PL_03_three.odum_std,'k')
text(PL_03_three.btl_time(7),-110,'C','FontSize', fontSize, 'FontName', fontName)

legend('Bottle NCP', 'PIGI NCP','Odum NCP','FontSize', fontSize, 'FontName', fontName)
% xlabel('Date','FontSize', fontSize, 'FontName', fontName)
% ylabel('NCP (mM/m^2/day)','FontSize', fontSize, 'FontName', fontName)
set(gca,'YTickLabel', [], 'FontSize', fontSize, 'FontName', fontName)
ylim([-125 75])
hold off

fig = gcf;
fig.Position = [100, 100, 1000, 500];

% Ensure high-quality rendering
set([ax1, ax2, ax3], 'FontSize', fontSize, 'FontName', fontName, 'LineWidth', lineWidth, 'Box', 'on')

% clear PL_03_three PL_02_three PL_01_three ax1 ax2 ax3 fig 
% print(fig, '-dpng', '-r600', 'NCP_comparison_3.png'); % PNG (for self)
% print(fig, 'NCP_comparison_3', '-depsc', '-r600'); % EPS (for publication)

%% Figure 6 - FFT in a 3x3 grid

% window = hann(L)'; % Hanning Window
% Y = fft(inter_O2 .* window); % Apply windowing before FFT


% PL_01_FFT
% Define base varables    
T = 10;                                   % Sampling period  (sec) 
Fs = 1/T;                                 % Sampling frequency (Hz)
L_PL01 = length(PL01.datetime);           % Length of signal
t = (0:L_PL01-1)*T;                       % Time vector

% need to fill missing data, using linerar interpolation
interp_NCP = fillmissing(PL01.NCP', 'linear');
interp_sal = fillmissing(PL01.qc_cal.sal', 'linear');
interp_solar = fillmissing(PL01.ts_cor.solar_rad', 'linear');
window_01 = hann(L_PL01)';
% Use the FFT
Y_NCP = fft(interp_NCP.* window_01);
power_spec_NCP_PL01 = abs(Y_NCP/L_PL01).^2;
Y_sal = fft(interp_sal.* window_01);
power_spec_sal_PL01 = abs(Y_sal/L_PL01).^2;
Y_solar = fft(interp_solar.* window_01);
power_spec_solar_PL01 = abs(Y_solar/L_PL01).^2;

f_PL_01 = Fs*(0:(L_PL01/2))/L_PL01;%frequency does not change between varables

% make 95% threshold
noise_NCP = randn(size(interp_NCP)) * mean(interp_NCP);
fft_noise_NCP = fft(noise_NCP.* window_01);
threhold_NCP_PL01 = prctile(abs(fft_noise_NCP/L_PL01).^2, 95);

noise_sal = randn(size(interp_sal)) * mean(interp_sal);
fft_noise_sal = fft(noise_sal.* window_01);
threhold_sal_PL01 = prctile(abs(fft_noise_sal/L_PL01).^2, 95);

noise_solar = randn(size(interp_solar)) * mean(interp_solar);
fft_noise_solar = fft(noise_solar.* window_01);
threhold_solar_PL01 = prctile(abs(fft_noise_solar/L_PL01).^2, 95);

clear fft_noise_solar fft_noise_sal fft_noise_NCP interp_solar interp_sal interp_NCP ...
    noise_solar noise_sal noise_NCP T t Y_solar Y_sal Y_NCP Fs window

%because there were cleanings, we need to adjust the time so that its
%continuious
continuios_time_PL02 = PL02.datetime(1):seconds(10):PL02.datetime(end-14);
%because datetime end is NaT, we moved 14 points back

% need to fill missing data, using linerar interpolation
full_NCP = fillmissing(PL02.NCP(1:end-14)', 'linear');
full_sal = fillmissing(PL02.qc_cal.sal(1:end-14)', 'linear');
full_solar = fillmissing(PL02.ts_cor.solar_rad(1:end-14)', 'linear');

time_pigi = PL02.datetime(1:end-14); %with out the last 14 NaT values

interp_NCP = interp1(time_pigi,full_NCP, continuios_time_PL02);
interp_sal = interp1(time_pigi,full_sal, continuios_time_PL02);
interp_solar = interp1(time_pigi,full_solar, continuios_time_PL02);

% PL_02_FFt
% Define base varables    
T = 10;                                   % Sampling period  (sec) 
Fs = 1/T;                                 % Sampling frequency (Hz)
L_PL02 = length(continuios_time_PL02);    % Length of signal
t = (0:L_PL02-1)*T;                       % Time vector
window_02 = hann(L_PL02)';
% Use the FFT
Y_NCP = fft(interp_NCP.*window_02);
power_spec_NCP_PL02 = abs(Y_NCP/L_PL02).^2;
Y_sal = fft(interp_sal.*window_02);
power_spec_sal_PL02 = abs(Y_sal/L_PL02).^2;
Y_solar = fft(interp_solar.*window_02);
power_spec_solar_PL02 = abs(Y_solar/L_PL02).^2;

f_PL_02 = Fs*(0:(L_PL02/2))/L_PL02;%frequency does not change between varables

% make 95% threshold
noise_NCP = randn(size(interp_NCP)) * mean(interp_NCP);
fft_noise_NCP = fft(noise_NCP.*window_02);
threhold_NCP_PL02 = prctile(abs(fft_noise_NCP/L_PL02).^2, 95);

noise_sal = randn(size(interp_sal)) * mean(interp_sal);
fft_noise_sal = fft(noise_sal.*window_02);
threhold_sal_PL02 = prctile(abs(fft_noise_sal/L_PL02).^2, 95);

noise_solar = randn(size(interp_solar)) * mean(interp_solar);
fft_noise_solar = fft(noise_solar.*window_02);
threhold_solar_PL02 = prctile(abs(fft_noise_solar/L_PL02).^2, 95);

clear fft_noise_solar fft_noise_sal fft_noise_NCP interp_solar interp_sal interp_NCP ...
    noise_solar noise_sal noise_NCP T t Y_solar Y_sal Y_NCP Fs full_NCP full_sal full_solar ...
    time_pigi continuios_time_PL02 window_02

% PL_03_FFt
% Define base varables    
T = 10;                                   % Sampling period  (sec) 
Fs = 1/T;                                 % Sampling frequency (Hz)
L_PL03 = length(PL03.datetime);           % Length of signal
t = (0:L_PL03-1)*T;                       % Time vector
window_03 = hann(L_PL03)';
% need to fill missing data, using linerar interpolation
interp_NCP = fillmissing(PL03.NCP', 'linear');
interp_sal = fillmissing(PL03.qc_cal.sal', 'linear');
interp_solar = fillmissing(PL03.ts_cor.solar_rad', 'linear');

% Use the FFT
Y_NCP = fft(interp_NCP.*window_03);
power_spec_NCP_PL03 = abs(Y_NCP/L_PL03).^2;
Y_sal = fft(interp_sal.*window_03);
power_spec_sal_PL03 = abs(Y_sal/L_PL03).^2;
Y_solar = fft(interp_solar.*window_03);
power_spec_solar_PL03 = abs(Y_solar/L_PL03).^2;

f_PL_03 = Fs*(0:(L_PL03/2))/L_PL03;%frequency does not change between varables

% make 95% threshold
noise_NCP = randn(size(interp_NCP)) * mean(interp_NCP);
fft_noise_NCP = fft(noise_NCP.*window_03);
threhold_NCP_PL03 = prctile(abs(fft_noise_NCP/L_PL03).^2, 95);

noise_sal = randn(size(interp_sal)) * mean(interp_sal);
fft_noise_sal = fft(noise_sal.*window_03);
threhold_sal_PL03 = prctile(abs(fft_noise_sal/L_PL03).^2, 95);

noise_solar = randn(size(interp_solar)) * mean(interp_solar);
fft_noise_solar = fft(noise_solar.*window_03);
threhold_solar_PL03 = prctile(abs(fft_noise_solar/L_PL03).^2, 95);

clear fft_noise_solar fft_noise_sal fft_noise_NCP interp_solar interp_sal interp_NCP ...
    noise_solar noise_sal noise_NCP T t Y_solar Y_sal Y_NCP Fs window_03

% plot figure

figure;
tiledlayout(3,3, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile;
plot(f_PL_01(1:floor(L_PL01/2)), power_spec_NCP_PL01(1:floor(L_PL01/2)), 'b', 'LineWidth', 1.5); % Power spectrum
hold on;
yline(threhold_NCP_PL01, 'r--', 'LineWidth', 1.5); % Significance threshold
text(4.45815*10^-5, 50, '6.23 hours', 'Rotation',90, 'FontSize', fontSize, 'FontName', fontName)
text(2.22907*10^-5, 160, '12.46 hours', 'Rotation',90, 'FontSize', fontSize, 'FontName', fontName)
text(1.11454*10^-5, 170, '24.92 hours', 'Rotation',90, 'FontSize', fontSize, 'FontName', fontName)
set(gca,'XTickLabel', [])
ylabel('Power - NCP')
text(0.9*10^-4,260,'A','FontSize', fontSize, 'FontName', fontName)
title('June')
hold off
xlim([0 1*10^-4])
ylim([0 300])

ax2 = nexttile;
plot(f_PL_02(1:floor(L_PL02/2)), power_spec_NCP_PL02(1:floor(L_PL02/2)), 'b', 'LineWidth', 1.5); % Power spectrum
hold on;
yline(threhold_NCP_PL02, 'r--', 'LineWidth', 1.5); % Significance threshold
text(4.51847*10^-5, 15, '6.14 hours', 'Rotation',90, 'FontSize', fontSize, 'FontName', fontName)
text(3.36258*10^-5, 20, '8.26 hours', 'Rotation',90, 'FontSize', fontSize, 'FontName', fontName)
text(2.20669*10^-5, 65, '12.58 hours', 'Rotation',90, 'FontSize', fontSize, 'FontName', fontName)
text(1.15589*10^-5, 50, '24.03 hours', 'Rotation',90, 'FontSize', fontSize, 'FontName', fontName)
set(gca,'XTickLabel', [])
text(0.9*10^-4,100,'B','FontSize', fontSize, 'FontName', fontName)
title('September')
hold off
xlim([0 1*10^-4])
ylim([0 120])

ax3 = nexttile;
plot(f_PL_03(1:floor(L_PL03/2)), power_spec_NCP_PL03(1:floor(L_PL03/2)), 'b', 'LineWidth', 1.5); % Power spectrum
hold on;
yline(threhold_NCP_PL03, 'r--', 'LineWidth', 1.5); % Significance threshold
text(4.50586*10^-5, 3, '6.16 hours', 'Rotation',90, 'FontSize', fontSize, 'FontName', fontName)
text(3.37939*10^-5, 2, '8.21 hours', 'Rotation',90, 'FontSize', fontSize, 'FontName', fontName)
text(2.25293*10^-5, 6, '12.32 hours', 'Rotation',90, 'FontSize', fontSize, 'FontName', fontName)
set(gca,'XTickLabel', [])
text(0.9*10^-4,8,'C','FontSize', fontSize, 'FontName', fontName)
title('December')
hold off
xlim([0 1*10^-4])
ylim([0 10])

ax4 = nexttile;
plot(f_PL_01(1:floor(L_PL01/2)), power_spec_sal_PL01(1:floor(L_PL01/2)), 'b', 'LineWidth', 1.5); % Power spectrum
hold on;
yline(threhold_sal_PL01, 'r--', 'LineWidth', 1.5); % Significance threshold
text(2.22907*10^-5, .18, '12.46 hours', 'FontSize', fontSize, 'FontName', fontName)
set(gca,'XTickLabel', [])
ylabel('Power - Salinity')
text(0.9*10^-4,.25,'D','FontSize', fontSize, 'FontName', fontName)
hold off
xlim([0 1*10^-4])
ylim([0 .3])

ax5 = nexttile;
plot(f_PL_02(1:floor(L_PL02/2)), power_spec_sal_PL02(1:floor(L_PL02/2)), 'b', 'LineWidth', 1.5); % Power spectrum
hold on;
yline(threhold_sal_PL02, 'r--', 'LineWidth', 1.5); % Significance threshold
text(2.20669*10^-5, .22, '12.58 hours', 'FontSize', fontSize, 'FontName', fontName)
set(gca,'XTickLabel', [])
text(0.9*10^-4,.25,'E','FontSize', fontSize, 'FontName', fontName)
hold off
xlim([0 1*10^-4])
ylim([0 .3])

ax6 = nexttile;
plot(f_PL_03(1:floor(L_PL03/2)), power_spec_sal_PL03(1:floor(L_PL03/2)), 'b', 'LineWidth', 1.5); % Power spectrum
hold on;
yline(threhold_sal_PL03, 'r--', 'LineWidth', 1.5); % Significance threshold
text(2.25293*10^-5, .18, '12.32 hours', 'FontSize', fontSize, 'FontName', fontName)
set(gca,'XTickLabel', [])
text(0.9*10^-4,.25,'F','FontSize', fontSize, 'FontName', fontName)
hold off
xlim([0 1*10^-4])
ylim([0 .3])

ax7 = nexttile;
plot(f_PL_01(1:floor(L_PL01/2)), power_spec_solar_PL01(1:floor(L_PL01/2)), 'b', 'LineWidth', 1.5); % Power spectrum
hold on;
yline(threhold_solar_PL01, 'r--', 'LineWidth', 1.5); % Significance threshold
text(1.11454*10^-5, 10000, '24.92 hours', 'FontSize', fontSize, 'FontName', fontName)
ylabel('Power - Solar Radation')
xlabel('Frequency (Hz)')
text(0.9*10^-4,18000,'G','FontSize', fontSize, 'FontName', fontName)
hold off
xlim([0 1*10^-4])
ylim([0 20000])

ax8 = nexttile;
plot(f_PL_02(1:floor(L_PL02/2)), power_spec_solar_PL02(1:floor(L_PL02/2)), 'b', 'LineWidth', 1.5); % Power spectrum
hold on;
yline(threhold_solar_PL02, 'r--', 'LineWidth', 1.5); % Significance threshold
text(2.31177*10^-5, 2500, '12.01 hours', 'FontSize', fontSize, 'FontName', fontName)
text(1.15589*10^-5, 9000, '24.03 hours', 'FontSize', fontSize, 'FontName', fontName)
xlabel('Frequency (Hz)')
text(0.9*10^-4,16000,'H','FontSize', fontSize, 'FontName', fontName)
hold off
xlim([0 1*10^-4])
ylim([0 18000])

ax9 = nexttile;
plot(f_PL_03(1:floor(L_PL03/2)), power_spec_solar_PL03(1:floor(L_PL03/2)), 'b', 'LineWidth', 1.5); % Power spectrum
hold on;
yline(threhold_solar_PL03, 'r--', 'LineWidth', 1.5); % Significance threshold
text(2.25293*10^-5, 600, '12.32 hours', 'FontSize', fontSize, 'FontName', fontName)
text(1.12646*10^-5, 1700, '24.65 hours', 'FontSize', fontSize, 'FontName', fontName)
xlabel('Frequency (Hz)')
text(0.9*10^-4,2700,'I','FontSize', fontSize, 'FontName', fontName)
hold off
xlim([0 1*10^-4])
ylim([0 3000])


fig = gcf;
fig.Position = [100, 100, 1000, 800];
% Ensure high-quality rendering
set([ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9], 'FontSize', fontSize, 'FontName', fontName, 'LineWidth', lineWidth, 'Box', 'on')

% print(fig, '-dpng', '-r600', 'FFT_03.png'); % PNG (for self)
% print(fig, 'FFT_03', '-depsc', '-r600'); % EPS (for publication)

clearvars -except PL02 fontSize fontName lineWidth PL01 PL03


%% Suplment 1 - Fix N2 gas

%define the good region to do the regression with
x = [ones(size(PL02.raw.o2uM(300:52650))), PL02.raw.o2uM(300:52650)];
y= PL02.raw.tp(300:52650);
%regress needs x to have a column of ones
[b,bint,r,rint,stats] = regress(y,x);
mdl = fitlm(PL02.raw.o2uM(300:52650),PL02.raw.tp(300:52650));
% giving the regression a try 
rf_tdg = PL02.raw.o2uM.*b(2)+b(1);

figure
hold on
plot(mdl)
xlabel('Oxygen (mM)', 'FontSize', fontSize, 'FontName', fontName)
ylabel('Total disolved Gas Pressure (mbar)', 'FontSize', fontSize, 'FontName', fontName)
text(125,1010, sprintf('R^2 = %.3f\n y = 0.752x + 846.94', stats(1)), 'FontSize', fontSize, 'FontName', fontName)
title('')
fig = gcf;
fig.Position = [100, 100, 1000, 800];
set(gca, 'FontSize', fontSize, 'FontName', fontName, 'LineWidth', lineWidth, 'Box', 'on')
hold off
% print(fig, '-dpng', '-r600', 'oxygen_to_TDGP.png'); % PNG (for self)
% print(fig, 'oxygen_to_TDGP', '-depsc', '-r600'); % EPS (for publication)

%%
figure
plot(PL02.datetime,PL02.raw.tp)
hold on
plot(PL02.datetime,rf_tdg)
legend('raw', 'fixed')
ylabel('Total disolved Gas Pressure (mbar)', 'FontSize', fontSize, 'FontName', fontName)
fig = gcf;
fig.Position = [100, 100, 1000, 800];
set(gca, 'FontSize', fontSize, 'FontName', fontName, 'LineWidth', lineWidth, 'Box', 'on')
hold off
% print(fig, '-dpng', '-r600', 'TDGP_fix.png'); % PNG (for self)
% print(fig, 'TDGP_fix', '-depsc', '-r600'); % EPS (for publication)

% clearvars -except PL02 fontSize fontName lineWidth PL01 PL03





