close all
clear all
clc

%Ben Lowin
%March 12th 2025
% estimating residuals between odums and various PIGI time windows.

%% Load data

load PL01_NCP_03.mat
PL01 = pigi_dat; % June data

load PL02_NCP_03.mat
PL02 = pigi_dat; % September data

load PL03_NCP_03.mat
PL03 = pigi_dat; % December data
clear pigi_dat

load PL03_three_plot.mat % June Odum
load PL02_three_plot.mat % September Odum
load PL01_three_plot.mat % December Odum

%% Base variables

time_around = 1:1:72; % Time window for averaging

% Initialize storage for residuals
r_time_June = nan(1,length(time_around));
r_time_Sept = nan(1,length(time_around));
r_time_Dec = nan(1,length(time_around));

% Define function to calculate residuals
calc_residuals = @(PIGI_data, odum_data) deal(...
    nan(1,length(time_around)), ...
    odum_data.odum_time, ...
    odum_data.odum_ncp, ...
    PIGI_data.datetime, ...
    PIGI_data.NCP ./ mean(PIGI_data.tide_depth,'omitnan'));

%% Compute residuals for each dataset

datasets = {PL01, PL02, PL03};
odum_datasets = {PL_01_three, PL_02_three, PL_03_three};
residual_results = {r_time_June, r_time_Sept, r_time_Dec};

for i = 1:3
    [residual_results{i}, odum_time, odum_NCP, PIGI_time, PIGI_NCP] = calc_residuals(datasets{i}, odum_datasets{i});
    
    for inn = 1:length(time_around)
        PIGI_mean = nan(length(odum_time),1);

        for in = 1:length(odum_time)
            top = odum_time(in) + hours(time_around(inn));
            bot = odum_time(in) - hours(time_around(inn));
            time_index = find(PIGI_time > bot & PIGI_time < top);
            PIGI_mean(in) = mean(PIGI_NCP(time_index), 'omitnan');
        end

        residual_results{i}(inn) = mean(PIGI_mean - odum_NCP', 'omitnan');
    end
end

% Assign computed residuals to specific variables
r_time_June = residual_results{1};
r_time_Sept = residual_results{2};
r_time_Dec = residual_results{3};

% Save results
save('Residuals_June.mat', 'r_time_June');
save('Residuals_Sept.mat', 'r_time_Sept');
save('Residuals_Dec.mat', 'r_time_Dec');

%% Plot results
r_mean = (r_time_Dec+r_time_Sept+r_time_June)./3;

figure
plot(time_around, r_time_June, '-r', 'LineWidth', 1.5); hold on;
plot(time_around, r_time_Sept, '-g', 'LineWidth', 1.5);
plot(time_around, r_time_Dec, '-b', 'LineWidth', 1.5);
plot(time_around,r_mean,'-k', 'LineWidth', 1.5)
hold off;

ylabel('NCP Residuals (mM/m^2/day)')
xlabel('Number of hours over which average is taken (h)')
set(gca, 'YDir', 'reverse')
legend({'June', 'September', 'December','Mean'}, 'Location', 'best')

fig = gcf;
fig.Position = [100, 100, 1000, 800];


% print(fig, '-dpng', '-r600', 'residual_pigi_odum_01.png'); % PNG (for self)
% print(fig, 'residual_pigi_odum_01', '-depsc', '-r600'); % EPS (for publication)


