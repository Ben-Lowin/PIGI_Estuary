close all;
clearvars;
clc;

% Ben Lowin
% March 18th, 2025
% Apply a phase lock to the NCP.

%% Load Data
load PL01_NCP_03.mat; PL01 = pigi_dat; % June data
load PL02_NCP_03.mat; PL02 = pigi_dat; % September data
load PL03_NCP_03.mat; PL03 = pigi_dat; % December data
clear pigi_dat;


%% Phase Lock Processing PL01
full_water = PL01.tide_depth;
interp_NCP = PL01.NCP./mean(full_water,'omitnan');
time_stamps = PL01.datetime; % Extract datetime vector

% Find local maxima and minima in tide depth
indexbot = find(islocalmax(full_water, 'MinProminence', 0.3, 'MinSeparation', 0.3, 'FlatSelection', 'first'));
indextop = find(islocalmin(full_water, 'MinProminence', 0.3, 'MinSeparation', 0.3, 'FlatSelection', 'first'));

a = floor(length(indextop)/2);


figure; hold on;
holder01 = cell(1, 5);
time_segments = cell(1, 5); % Store corresponding time vectors
inn=1;

for in = 1:5
    segment = interp_NCP(indextop(inn):indextop(inn+2));
    time_segment = time_stamps(indextop(inn):indextop(inn+2)); % Extract time
    
    plot(time_segment, segment);
    holder01{in} = segment;
    time_segments{in} = time_segment;
    inn=inn+2;
end

hold off;
xlabel('Time of Tidal Phase');
ylabel('NCP');
title('Tidal Phase Segments');
grid on;
datetick('x', 'HH:MM', 'keeplimits'); % Format x-axis to show hours and minutes

% Interpolation to Common Length
max_len = max(cellfun(@length, holder01));
interp_segments = NaN(length(holder01), max_len);

for i = 1:length(holder01)
    segment = holder01{i};
    time_segment = time_segments{i};
    
    % Normalize time to a common scale (e.g., percentage of full phase)
    x_old = linspace(0, 1, length(segment));
    x_new = linspace(0, 1, max_len);
    
    % Interpolate NCP values
    interp_segments(i, :) = interp1(x_old, segment, x_new, 'linear', 'extrap');
end

% Compute mean while ignoring NaNs
avg_segment = mean(interp_segments, 1, 'omitnan');

% Plot results
figure; hold on;
for i = 1:size(interp_segments, 1)
    plot(linspace(0, 1, max_len), interp_segments(i, :), 'Color', [0.7, 0.7, 0.7]); % Light gray
end
plot(linspace(0, 1, max_len), avg_segment, 'k', 'LineWidth', 2); % Bold black for mean

xlabel('Normalized Time of Tidal Phase');
ylabel('NCP');
title('NCP Segments and The Average');
grid on;
legend({'Base Segments', 'Average NCP'}, 'Location', 'Best');
hold off;
%% Phase Lock Processing PL02
full_water = PL02.tide_depth;
interp_NCP = PL02.NCP./mean(full_water,'omitnan');
time_stamps = PL02.datetime; % Extract datetime vector

% Find local maxima and minima in tide depth
indexbot = find(islocalmax(full_water, 'MinProminence', 0.3, 'MinSeparation', 0.3, 'FlatSelection', 'first'));
indextop = find(islocalmin(full_water, 'MinProminence', 0.3, 'MinSeparation', 0.3, 'FlatSelection', 'first'));

figure; hold on;
holder02 = cell(1, 10);
time_segments = cell(1, 10); % Store corresponding time vectors
inn=1;

for in = 1:10
    segment = interp_NCP(indextop(inn):indextop(inn+2));
    time_segment = time_stamps(indextop(inn):indextop(inn+2)); % Extract time
    
    plot(time_segment, segment);
    holder02{in} = segment;
    time_segments{in} = time_segment;
    inn=inn+2;
end

hold off;
xlabel('Time of Tidal Phase');
ylabel('NCP');
title('Tidal Phase Segments');
grid on;
datetick('x', 'HH:MM', 'keeplimits'); % Format x-axis to show hours and minutes

% Interpolation to Common Length
max_len = max(cellfun(@length, holder02));
interp_segments = NaN(length(holder02), max_len);

for i = 1:length(holder02)
    segment = holder02{i};
    time_segment = time_segments{i};
    
    % Normalize time to a common scale (e.g., percentage of full phase)
    x_old = linspace(0, 1, length(segment));
    x_new = linspace(0, 1, max_len);
    
    % Interpolate NCP values
    interp_segments(i, :) = interp1(x_old, segment, x_new, 'linear', 'extrap');
end

% Compute mean while ignoring NaNs
avg_segment = mean(interp_segments, 1, 'omitnan');

% Plot results
figure; hold on;
for i = 1:size(interp_segments, 1)
    plot(linspace(0, 1, max_len), interp_segments(i, :), 'Color', [0.7, 0.7, 0.7]); % Light gray
end
plot(linspace(0, 1, max_len), avg_segment, 'k', 'LineWidth', 2); % Bold black for mean

xlabel('Normalized Time of Tidal Phase');
ylabel('NCP');
title('NCP Segments and The Average');
grid on;
legend({'Base Segments', 'Average NCP'}, 'Location', 'Best');
hold off;

%% Phase Lock Processing PL03
full_water = PL03.tide_depth;
interp_NCP = PL03.NCP./mean(full_water,'omitnan');
time_stamps = PL03.datetime; % Extract datetime vector

% Find local maxima and minima in tide depth
indexbot = find(islocalmax(full_water, 'MinProminence', 0.3, 'MinSeparation', 0.3, 'FlatSelection', 'first'));
indextop = find(islocalmin(full_water, 'MinProminence', 0.3, 'MinSeparation', 0.3, 'FlatSelection', 'first'));

figure; hold on;
holder03 = cell(1, 5);
time_segments = cell(1, 5); % Store corresponding time vectors
inn=1;
for in = 1:5
    segment = interp_NCP(indextop(inn):indextop(inn+2));
    time_segment = time_stamps(indextop(inn):indextop(inn+2)); % Extract time
    
    plot(time_segment, segment);
    holder03{in} = segment;
    time_segments{in} = time_segment;
   inn=inn+2;
end

hold off;
xlabel('Time of Tidal Phase');
ylabel('NCP');
title('Tidal Phase Segments');
grid on;
datetick('x', 'HH:MM', 'keeplimits'); % Format x-axis to show hours and minutes

% Interpolation to Common Length
max_len = max(cellfun(@length, holder03));
interp_segments = NaN(length(holder03), max_len);

for i = 1:length(holder03)
    segment = holder03{i};
    time_segment = time_segments{i};
    
    % Normalize time to a common scale (e.g., percentage of full phase)
    x_old = linspace(0, 1, length(segment));
    x_new = linspace(0, 1, max_len);
    
    % Interpolate NCP values
    interp_segments(i, :) = interp1(x_old, segment, x_new, 'linear', 'extrap');
end

% Compute mean while ignoring NaNs
avg_segment = mean(interp_segments, 1, 'omitnan');

% Plot results
figure; hold on;
for i = 1:size(interp_segments, 1)
    plot(linspace(0, 1, max_len), interp_segments(i, :), 'Color', [0.7, 0.7, 0.7]); % Light gray
end
plot(linspace(0, 1, max_len), avg_segment, 'k', 'LineWidth', 2); % Bold black for mean

xlabel('Normalized Time of Tidal Phase');
ylabel('NCP');
title('NCP Segments and The Average');
grid on;
legend({'Base Segments', 'Average NCP'}, 'Location', 'Best');
hold off;
