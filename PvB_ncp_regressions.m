close all
clear all
clc

% Ben Lowin
% July 26th 2024
% Regressions for PIGI vs bottle data

%%
%pigi data
load PL01_NCP_01.mat
%bottle data + ncp pigi in m^3
load NCP_bottle_data.mat
%water height/tide
load Water_height.mat
%YSI/ancilary data
load YSI_data.mat
%Weather data
load Weather_data.mat
%hanna_data
load Hanna_data.mat

%%

hanna_datetime = Hanna.date+Hanna.time;

%% interp all data to be on the same time ticks

%NCP
pvb_regress.datetime = pigi_dat.datetime;
pvb_regress.NCP = bottle.pigi_NCP_m3;
%Hanna
pvb_regress.temp = interp1(hanna_datetime,Hanna.temp,pvb_regress.datetime);
pvb_regress.salinity = interp1(hanna_datetime,Hanna.salinity,pvb_regress.datetime);
pvb_regress.pH = interp1(hanna_datetime,Hanna.pH,pvb_regress.datetime);
pvb_regress.DO = interp1(hanna_datetime,Hanna.Disolved_oxygen,pvb_regress.datetime);
%weather
pvb_regress.wind_speed = interp1(weather.datetime,weather.wind_speed,pvb_regress.datetime);
pvb_regress.air_temp = interp1(weather.datetime,weather.air_temp,pvb_regress.datetime);
pvb_regress.preasure = interp1(weather.datetime,weather.preasure,pvb_regress.datetime);
pvb_regress.PAR = interp1(weather.datetime,weather.PAR,pvb_regress.datetime);
%water_height
pvb_regress.water_depth = interp1(water_height.time,water_height.depth_real,pvb_regress.datetime);
pvb_regress.NCP_ones = [ones(size(pvb_regress.NCP,1),1),pvb_regress.NCP];
  
clear a bottle Hanna hanna_datetime pigi_dat water_height weather YSI

%% regression for all points
% % [B,BINT,R,RINT,STATS] = regress (pvb_regress.temp,pvb_regress.NCP_ones);
% 
% % Variables to loop through for regression
% variables = {'temp', 'salinity', 'pH', 'DO', 'wind_speed', 'air_temp', 'preasure', 'PAR', 'water_depth'};
% 
% % Initialize results storage
% results = struct();
% 
% % Loop through each variable and perform regression
% for i = 1:length(variables)
%     var_name = variables{i};
%     Y = pvb_regress.(var_name);
%     
%     % Perform regression
%     [B, BINT, R, RINT, STATS] = regress(Y, pvb_regress.NCP_ones);
%     
%     % Store results
%     results.(var_name).B = B;
%     results.(var_name).BINT = BINT;
%     results.(var_name).R = R;
%     results.(var_name).RINT = RINT;
%     results.(var_name).STATS = STATS;
% 
%     % Plot the regression line and data points
%     figure;
%     scatter(pvb_regress.NCP, Y, '.');
%     hold on;
%     % Regression line
%     x_fit = linspace(min(pvb_regress.NCP), max(pvb_regress.NCP), 100);
%     y_fit = B(1) + B(2)*x_fit;
%     plot(x_fit, y_fit, '-r', 'LineWidth', 2);
%     hold off;
%     title(sprintf('Regression of %s on NCP', var_name));
%     xlabel('NCP');
%     ylabel(var_name);
%     legend('Data points', 'Regression line');
% end
% % STATS containing, in the following order, the R-square statistic,
% % the F statistic and p value for the full model, and an estimate of the error variance.
% 
% clear B BINT i R RINT STATS var_name Y

%%

load NCP_bottle_data.mat


min(pvb_regress.PAR)

for in =1:10

time_out=bottle.sample_dates(in)+hours(24);

PAR_index = find(pvb_regress.datetime>bottle.sample_dates(in) & pvb_regress.datetime< time_out);

mean_PAR(in) = mean(pvb_regress.PAR(PAR_index));

end

plot (mean_PAR,bottle.all_NCP,'.')
xlabel('PAR')
ylabel('NCP (umol O_2/day)')

