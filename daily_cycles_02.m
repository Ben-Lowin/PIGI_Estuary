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
fontSize = 12;
fontName = 'Times New Roman';
lineWidth = 1.2;


%% Process Data
[interp_June, avg_June, max_len_June] = process_NCP(PL01.datetime,...
    PL01.NCP./mean(PL01.tide_depth,'omitnan'));
[interp_Sept, avg_Sept, max_len_Sept] = process_NCP(PL02.datetime,...
    PL02.NCP./mean(PL02.tide_depth,'omitnan'));
[interp_Dec, avg_Dec, max_len_Dec] = process_NCP(PL03.datetime,...
    PL03.NCP./mean(PL03.tide_depth,'omitnan'));

%% Create 2x3 Panel Figure
figure;

% June - Segments
subplot(3,1,1); hold on;
for i = 1:size(interp_June, 1)
    plot(linspace(0, 1, max_len_June), interp_June(i, :), 'Color', [0.7, 0.7, 0.7]); % Light gray
end
ylabel('NCP (mmol O_2 /m^2/day)','FontSize',fontSize, 'FontName', fontName);
xticks([0, 0.5, 1])
xticklabels(['6am'; '6pm'; '6am'])
title('June','FontSize',fontSize, 'FontName', fontName);
plot(linspace(0, 1, max_len_June), avg_June, 'k', 'LineWidth', 2);
text(0.05,-25,'A','FontSize',fontSize, 'FontName', fontName)
ylim([-250 0])
hold off;

% September -  Segments
subplot(3,1,2); hold on;
for i = 1:size(interp_Sept, 1)
    plot(linspace(0, 1, max_len_Sept), interp_Sept(i, :), 'Color', [0.7, 0.7, 0.7]); % Light gray
end
ylabel('NCP (mmol O_2 /m^2/day)','FontSize',fontSize, 'FontName', fontName);
title('September','FontSize',fontSize, 'FontName', fontName);
xticks([0, 0.5, 1])
xticklabels(['6am'; '6pm'; '6am'])
plot(linspace(0, 1, max_len_Sept), avg_Sept, 'k', 'LineWidth', 2);
text(0.05,-25,'B','FontSize',fontSize, 'FontName', fontName)
ylim([-250 0])
hold off;

% December - Individual Segments
subplot(3,1,3); hold on;
for i = 1:size(interp_Dec, 1)
    plot(linspace(0, 1, max_len_Dec), interp_Dec(i, :), 'Color', [0.7, 0.7, 0.7]); % Light gray
end
xlabel('Time of Day (EST)','FontSize',fontSize, 'FontName', fontName);
ylabel('NCP (mmol O_2 /m^2/day)','FontSize',fontSize, 'FontName', fontName);
title('December','FontSize',fontSize, 'FontName', fontName);
xticks([0, 0.5, 1])
xticklabels(['6am'; '6pm'; '6am'])
plot(linspace(0, 1, max_len_Dec), avg_Dec, 'k', 'LineWidth', 2);
text(0.05,-15,'C','FontSize',fontSize, 'FontName', fontName)
ylim([-150 0])
hold off;

fig = gcf;
fig.Position = [100, 100, 600, 800];


% print(fig, '-dpng', '-r600', 'daily_cycle_01.png'); % PNG (for self)
% print(fig, 'daily_cycle_01', '-depsc', '-r600'); % EPS (for publication)


%% Define Function for Processing and Plotting
function [interp_segments, avg_segment, max_len] = process_NCP(time, NCP)
    day_point = hour(time);
    day_index= find(day_point==0);
    inin = find(day_index(2:end)-day_index(1:end-1)>2);
    ticksi = day_index(inin);

    hold_data = {};
    time_segments = {};

    for in = 1:length(inin)-1
        segment = NCP(ticksi(in):ticksi(in+1));
        time_segment = time(ticksi(in):ticksi(in+1));

        hold_data{in} = segment;
        time_segments{in} = time_segment;
    end

    % Interpolation to Common Length
    max_len = max(cellfun(@length, hold_data));
    interp_segments = NaN(length(hold_data), max_len);

    for i = 1:length(hold_data)
        segment = hold_data{i};
        x_old = linspace(0, 1, length(segment));
        x_new = linspace(0, 1, max_len);
        interp_segments(i, :) = interp1(x_old, segment, x_new, 'linear', 'extrap');
    end

    % Compute mean while ignoring NaNs
    avg_segment = mean(interp_segments, 1, 'omitnan');
end
