close all;
clearvars;
clc;

% Ben Lowin
% March 18th, 2025
% Partition data into 24-hour segments from 6 AM to 6 AM.

%% Load Data
load PL01_NCP_03.mat; PL01 = pigi_dat; % June data
load PL02_NCP_03.mat; PL02 = pigi_dat; % September data
load PL03_NCP_03.mat; PL03 = pigi_dat; % December data
clear pigi_dat;


%%
time = PL01.datetime;
NCP = PL01.NCP;

day_point = hour(time);
day_index= find(day_point==0);
inin = find(day_index(2:end)-day_index(1:end-1)>2);
ticksi = day_index(inin);

figure
hold on

for in = 1:length(inin)-1

    segment = NCP(ticksi(in):ticksi(in+1));
    time_segment = time(ticksi(in):ticksi(in+1));

    plot(time_segment, segment);
    holder01{in} = segment;
    time_segments{in} = time_segment;

 end
hold off

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

%%
time = PL02.datetime;
NCP = PL02.NCP;

day_point = hour(time);
day_index= find(day_point==0);
inin = find(day_index(2:end)-day_index(1:end-1)>2);
ticksi = day_index(inin);

figure
hold on

for in = 1:length(inin)-1

    segment = NCP(ticksi(in):ticksi(in+1));
    time_segment = time(ticksi(in):ticksi(in+1));

    plot(time_segment, segment);
    holder02{in} = segment;
    time_segments{in} = time_segment;

 end
hold off

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

%%

time = PL03.datetime;
NCP = PL03.NCP;

day_point = hour(time);
day_index= find(day_point==0);
inin = find(day_index(2:end)-day_index(1:end-1)>2);
ticksi = day_index(inin);

figure
hold on

for in = 1:length(inin)-1

    segment = NCP(ticksi(in):ticksi(in+1));
    time_segment = time(ticksi(in):ticksi(in+1));

    plot(time_segment, segment);
    holder03{in} = segment;
    time_segments{in} = time_segment;

 end
hold off

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

