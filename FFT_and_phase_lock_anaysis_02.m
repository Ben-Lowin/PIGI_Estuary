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

%% Formatting Variables
fontSize = 12;
fontName = 'Times New Roman';
lineWidth = 1.2;

%% Process Data
[interp_June, avg_June, max_len_June] = process_NCP(PL01);
[interp_Sept, avg_Sept, max_len_Sept] = process_NCP(PL02);
[interp_Dec, avg_Dec, max_len_Dec] = process_NCP(PL03);

%% Create 3-Panel Figure (Tidal Cycle)
figure;

% June - Tidal Phase Segments
subplot(3,1,1); hold on;
for i = 1:size(interp_June, 1)
    plot(linspace(0, 1, max_len_June), interp_June(i, :), 'Color', [0.7, 0.7, 0.7]); % Light gray
end
ylabel('NCP (mmol O_2 /m^2/day)', 'FontSize', fontSize, 'FontName', fontName);
% ylabel('water depth (m)', 'FontSize', fontSize, 'FontName', fontName);
% ylabel('Oxygen (mM)', 'FontSize', fontSize, 'FontName', fontName);
xticks([0, 0.5, 1]);
xticklabels({'Low Tide','High Tide', 'Low Tide'});
title('June', 'FontSize', fontSize, 'FontName', fontName);
plot(linspace(0, 1, max_len_June), avg_June, 'k', 'LineWidth', 2);
text(0.05, -25, 'A', 'FontSize', fontSize, 'FontName', fontName);
% ylim([-150 0])
hold off;

% September - Tidal Phase Segments
subplot(3,1,2); hold on;
for i = 1:size(interp_Sept, 1)
    plot(linspace(0, 1, max_len_Sept), interp_Sept(i, :), 'Color', [0.7, 0.7, 0.7]); % Light gray
end
ylabel('NCP (mmol O_2 /m^2/day)', 'FontSize', fontSize, 'FontName', fontName);
xticks([0, 0.5, 1]);
xticklabels({'Low Tide','High Tide', 'Low Tide'});
title('September', 'FontSize', fontSize, 'FontName', fontName);
plot(linspace(0, 1, max_len_Sept), avg_Sept, 'k', 'LineWidth', 2);
text(0.05, -25, 'B', 'FontSize', fontSize, 'FontName', fontName);
% ylim([-150 0])
hold off;

% December - Tidal Phase Segments
subplot(3,1,3); hold on;
for i = 1:size(interp_Dec, 1)
    plot(linspace(0, 1, max_len_Dec), interp_Dec(i, :), 'Color', [0.7, 0.7, 0.7]); % Light gray
end
xlabel('Tidal Phase', 'FontSize', fontSize, 'FontName', fontName);
ylabel('NCP (mmol O_2 /m^2/day)', 'FontSize', fontSize, 'FontName', fontName);
xticks([0, 0.5, 1]);
xticklabels({'Low Tide','High Tide', 'Low Tide'});
title('December', 'FontSize', fontSize, 'FontName', fontName);
plot(linspace(0, 1, max_len_Dec), avg_Dec, 'k', 'LineWidth', 2);
text(0.05, -15, 'C', 'FontSize', fontSize, 'FontName', fontName);
% ylim([-100 0])
hold off;

fig = gcf;
fig.Position = [100, 100, 600, 800];

% Save as PNG and EPS
% print(fig, '-dpng', '-r600', 'tidal_cycle_01.png'); % PNG (for self)
% print(fig, 'tidal_cycle_01', '-depsc', '-r600'); % EPS (for publication)

%% Function for Processing Tidal Phase NCP
function [interp_segments, avg_segment, max_len] = process_NCP(data)
    full_water = data.tide_depth;
%     interp_NCP = full_water;
    interp_NCP = data.NCP ./ mean(full_water, 'omitnan');
%     interp_NCP = data.ts_cor.o2uM;
    time_stamps = data.datetime;

    % Find local maxima and minima in tide depth
    indextop = find(islocalmin(full_water, 'MinProminence', 0.3, 'MinSeparation', 0.3, 'FlatSelection', 'first'));

    % Segment extraction
    num_segments = length(indextop);
    holder = cell(1, num_segments);
    time_segments = cell(1, num_segments);
    inn = 1;

    for in = 1:num_segments-1
        segment = interp_NCP(indextop(inn):indextop(inn+1));
        time_segment = time_stamps(indextop(inn):indextop(inn+1));
        holder{in} = segment;
        time_segments{in} = time_segment;
        inn = inn + 1;
    end

    % Interpolation to Common Length
    max_len = max(cellfun(@length, holder));
    interp_segments = NaN(length(holder), max_len);

    for i = 1:length(holder)-1
        segment = holder{i};
        x_old = linspace(0, 1, length(segment));
        x_new = linspace(0, 1, max_len);
        interp_segments(i, :) = interp1(x_old, segment, x_new, 'linear', 'extrap');
    end

    % Compute mean while ignoring NaNs
    avg_segment = mean(interp_segments, 1, 'omitnan');
end
