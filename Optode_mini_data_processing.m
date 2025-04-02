close all
clear all
clc

% Ben Lowin
% June 28th 2024
% calculate the NCP, PP and respiration for each incubation
% determin the YSI means for each station
% save out to a single insitu data file

%% load in the incubation data
load Optode_mini_data.mat


%% Loop to calculate the mean, median and std of each sample

% Get the list of main fields in Optode_mini
dateFields = fieldnames(Optode_mini);

k=0;

% Loop through each date field
for i = 1:length(dateFields)
    dateField = dateFields{i};
    
    % Skip the 'units' field
    if strcmp(dateField, 'units')
        continue;
    end
    
    % Get the list of subfields for the current date field
    subFields = fieldnames(Optode_mini.(dateField));
    
    % Loop through each subfield
    for j = 1:length(subFields)
        subField = subFields{j};
        
        % Extract the oxygen data
        oxygenData = Optode_mini.(dateField).(subField).oxygen;
        k=k+1;
        % Calculate mean, median, and standard deviation
        sample_label (k) = {subField};
        meanOxygen(k) = mean(oxygenData);
        medianOxygen(k) = median(oxygenData);
        stdOxygen(k) = std(oxygenData);
    end
end

clear i j k dateField dateFields oxygenData subField subFields
%%  for each sample calcualte the NCP, repiration

% Initialize variables to store results
low_positive_means = [];
low_positive_stds = [];
low_negative_means = [];
low_negative_stds = [];
high_positive_means = [];
high_positive_stds = [];
high_negative_means = [];
high_negative_stds = [];

for startIndex = 1:24:120
    % Define the current range of indices
    i = startIndex:2:(startIndex+22);
    ii = (startIndex+1):2:(startIndex+23);

    % Calculate the change in oxygen
    change_in_oxygen = meanOxygen(i) - meanOxygen(ii);

    % Identify low tide samples
    low_tide_i = contains(sample_label(i), 'low');

    % Low tide oxygen changes
    low_tide_oxygen = change_in_oxygen(low_tide_i);
    low_positive = low_tide_oxygen(low_tide_oxygen > 0);
    low_negative = low_tide_oxygen(low_tide_oxygen < 0);

    % Calculate low tide statistics
    low_positive_mean = mean(low_positive);
    low_positive_std = std(low_positive);
    low_negative_mean = mean(low_negative);
    low_negative_std = std(low_negative);

    % Store low tide results
    low_positive_means = [low_positive_means, low_positive_mean]; %#ok<AGROW>
    low_positive_stds = [low_positive_stds, low_positive_std]; %#ok<AGROW>
    low_negative_means = [low_negative_means, low_negative_mean]; %#ok<AGROW>
    low_negative_stds = [low_negative_stds, low_negative_std]; %#ok<AGROW>

    % High tide oxygen changes
    high_tide_oxygen = change_in_oxygen(~low_tide_i);
    high_positive = high_tide_oxygen(high_tide_oxygen > 0);
    high_negative = high_tide_oxygen(high_tide_oxygen < 0);

    % Calculate high tide statistics
    high_positive_mean = mean(high_positive);
    high_positive_std = std(high_positive);
    high_negative_mean = mean(high_negative);
    high_negative_std = std(high_negative);

    % Store high tide results
    high_positive_means = [high_positive_means, high_positive_mean]; %#ok<AGROW>
    high_positive_stds = [high_positive_stds, high_positive_std]; %#ok<AGROW>
    high_negative_means = [high_negative_means, high_negative_mean]; %#ok<AGROW>
    high_negative_stds = [high_negative_stds, high_negative_std]; %#ok<AGROW>
end

clear i ii change_in_oxygen low_tide_oxygen low_positive low_negative ...
    low_positive_std low_positive_mean low_negative_std low_negative_mean ...
    high_tide_oxygen high_positive high_negative high_positive_mean high_positive_std ...
    high_negative_mean high_negative_std low_tide_i startIndex 

%% relable

bottle.HT_respiration_mean = high_negative_means;
bottle.HT_respiration_std = high_negative_stds;
bottle.HT_NCP_mean = high_positive_means;
bottle.HT_NCP_std = high_positive_stds;

bottle.LT_respiration_mean = low_negative_means;
bottle.LT_respiration_std = low_negative_stds;
bottle.LT_NCP_mean = low_positive_means;
bottle.LT_NCP_std = low_positive_stds;

clear high_negative_means high_negative_stds high_positive_means high_positive_stds ...
    low_negative_means low_positive_means low_positive_stds low_negative_stds
%% 
% NCP = photosynetesis - respiration
% change in light bottle = unknown - change in dark bottle
% delta light + delta dark = X

bottle.HT_photo_mean = bottle.HT_NCP_mean - bottle.HT_respiration_mean;
bottle.HT_photo_std = bottle.HT_NCP_std - bottle.HT_respiration_std;
bottle.LT_photo_mean = bottle.LT_NCP_mean - bottle.LT_respiration_mean;
bottle.LT_photo_std = bottle.LT_NCP_std - bottle.LT_respiration_std;


%%

figure 
errorbar(bottle.HT_respiration_mean,bottle.HT_respiration_std,'r*')
hold on
errorbar(bottle.HT_photo_mean,bottle.HT_photo_std,'g*')
errorbar(bottle.HT_NCP_mean,bottle.HT_NCP_std,'b*')
xlim([0 6])
xticklabels({'', '06/17', '06/18', '06/19', '06/20', '06/21', ''});
ylabel('mmol O_2/m^3/day')
title('High Tide incubations')
legend('repiration','photosynthesis','NCP')
hold off

figure 
errorbar(bottle.LT_respiration_mean,bottle.LT_respiration_std,'r*')
hold on
errorbar(bottle.LT_photo_mean,bottle.LT_photo_std,'g*')
errorbar(bottle.LT_NCP_mean,bottle.LT_NCP_std,'b*')
xlim([0 6])
xticklabels({'', '06/17', '06/18', '06/19', '06/20', '06/21', ''});
ylabel('mmol O_2/m^3/day')
title('Low Tide incubations')
legend('repiration','photosynthesis','NCP')
hold off

figure 
errorbar(bottle.HT_respiration_mean,bottle.HT_respiration_std,'r*')
hold on
errorbar(bottle.LT_respiration_mean,bottle.LT_respiration_std,'b*')
xlim([0 6])
xticklabels({'', '06/17', '06/18', '06/19', '06/20', '06/21', ''});
ylabel('\mumol O_2/L/day')
title('Respiration')
legend('High tide', 'Low tide')
hold off

figure 
errorbar(bottle.HT_photo_mean,bottle.HT_photo_std,'r*')
hold on
errorbar(bottle.LT_photo_mean,bottle.LT_photo_std,'b*')
xlim([0 6])
xticklabels({'', '06/17', '06/18', '06/19', '06/20', '06/21', ''});
ylabel('\mumol O_2/L/day')
title('Photosynthesis')
legend('High tide', 'Low tide')
hold off

figure 
errorbar(bottle.HT_NCP_mean,bottle.HT_NCP_std,'r*')
hold on
errorbar(bottle.LT_NCP_mean,bottle.LT_NCP_std,'b*')
xlim([0 6])
xticklabels({'', '06/17', '06/18', '06/19', '06/20', '06/21', ''});
ylabel('\mumol O_2/L/day')
title('Net Comunity Production')
legend('High tide', 'Low tide')
hold off

%% load NCP for comparison
load PL01_NCP_01.mat

figure
plot(pigi_dat.datetime,pigi_dat.NCP)
xlabel('Date')
ylabel('NCP (umol O_2/m^2/day)')

figure
plot(pigi_dat.datetime,pigi_dat.n2_calc.do2n2)
xlabel('Date')
ylabel('delta O_2/N_2 (uM O_2/day)')


%% put date time in bottles so they can be plotted
bottle.sample_dates(1) = datetime(2024,6,17,5,30,00,00);
bottle.sample_dates(2) = datetime(2024,6,17,11,53,00,00);
bottle.sample_dates(3) = datetime(2024,6,18,6,20,00,00);
bottle.sample_dates(4) = datetime(2024,6,18,12,50,00,00);
bottle.sample_dates(5) = datetime(2024,6,19,7,10,00,00);
bottle.sample_dates(6) = datetime(2024,6,19,13,23,00,00);
bottle.sample_dates(7) = datetime(2024,6,20,7,55,00,00);
bottle.sample_dates(8) = datetime(2024,6,20,14,20,00,00);
bottle.sample_dates(9) = datetime(2024,6,21,8,45,00,00);
bottle.sample_dates(10) = datetime(2024,6,21,15,03,00,00);

count = 0;
for in = 1:2:10
    count = count+1;
bottle.all_NCP(in)  = bottle.HT_NCP_mean(count);
bottle.all_NCP(in+1)  = bottle.LT_NCP_mean(count);
bottle.all_NCP_std(in)  = bottle.HT_NCP_std(count);
bottle.all_NCP_std(in+1)  = bottle.LT_NCP_std(count);
end


%%

NCP_m3=pigi_dat.NCP./3.5;

i= 1;

for i = 1:10

range = find (pigi_dat.datetime < bottle.sample_dates(i) + hours(1) & ...
    pigi_dat.datetime > bottle.sample_dates(i) - hours(1));

NCP_avr(i)=mean(NCP_m3(range),"omitnan");
NCP_std(i)=std(NCP_m3(range),"omitnan");

clear range
end


figure
plot(pigi_dat.datetime,NCP_m3)
xlabel('Date')
ylabel('NCP (mmol O_2/m^3/day)')

%%

figure
plot(1:10, NCP_avr,'*')
errorbar(1:10,NCP_avr,NCP_std,NCP_std,'*')
hold on
errorbar(1:10,bottle.all_NCP,bottle.all_NCP_std,bottle.all_NCP_std,'*')
xticklabels({'M-H','M-L','T-H','T-L','W-H','W-L','Th-H','Th-L','F-H','F-L'})
xlim([0 11])
ylabel('NCP pigi (mM O_2/m^3/day)')
legend ('Pigi data','Bottle data')


difference_bottle_pigi = NCP_avr - bottle.all_NCP;

bottle.pigi_ncp_avr = NCP_avr;
bottle.pigi_ncp_std = NCP_std;
bottle.pigi_NCP_m3 = NCP_m3;
bottle.pigi_datetime = pigi_dat.datetime;

% save('NCP_bottle_data','bottle')


% Bottle_m3 = bottle.all_NCP
%%
load Water_height.mat
PL_water_depth = water_height.depth +1.9;

% Convert the data to a timetable
TT = timetable(pigi_dat.datetime, NCP_m3);

% Resample the data to hourly intervals using mean to smooth
TT_hourly = retime(TT, 'hourly', 'mean');

% Extract the datetime and smoothed NCP values if needed
datetime_hourly = TT_hourly.datetime;
NCP_hourly = TT_hourly.NCP_m3;




figure
ax(1) = subplot(2,1,1);
plot(pigi_dat.datetime,NCP_m3)
xlabel('Date')
ylabel('NCP pigi (mM O_2/m^3/day)')

ax(2) = subplot(2,1,2);
plot(water_height.time,PL_water_depth)
xlabel('Date')
ylabel('Depth of water (m)')

linkaxes(ax,'x')


%%
% Ensure both datetime vectors are in the same format
pigi_dat.datetime = datetime(pigi_dat.datetime);
water_height.time = datetime(water_height.time);

% Remove non-finite values
finiteIndices = isfinite(pigi_dat.datetime) & isfinite(NCP_m3);
clean_datetime = pigi_dat.datetime(finiteIndices);
clean_NCP_m3 = pigi_dat.qc_cal.o2uM(finiteIndices);

% Perform the interpolation
NCP_water_interp = interp1(clean_datetime, clean_NCP_m3, water_height.time);

figure
plot(water_height.depth_real,NCP_water_interp,'.')

regress(NCP_water_interp,water_height.depth_real)

%%

figure
yyaxis left
plot(pigi_dat.datetime,NCP_m3)
xlabel('Date')
ylabel('NCP pigi (mM O_2/m^3/day)')
ylim([-140 80])

yyaxis right
errorbar(bottle.sample_dates,bottle.all_NCP,bottle.all_NCP_std,bottle.all_NCP_std,'*')
xlabel('Date')
ylabel('NCP bottle (mM O_2/m^3/day)')
ylim([-140 80])


% %% Plot Hanna data
% load Hanna_data.mat
% Hanna.date_time= Hanna.date+ Hanna.time; %make datatime
% 
% figure
% plot(Hanna.date_time,Hanna.temp)
% xlabel('Date')
% ylabel('Temperature (c)')
% 
% figure
% plot(Hanna.date_time,Hanna.salinity)
% xlabel('Date')
% ylabel('Salinity (PSU)')
% 
% figure
% plot(Hanna.date_time,Hanna.pH)
% xlabel('Date')
% ylabel('pH')
% 
% figure
% plot(Hanna.date_time,Hanna.Disolved_oxygen)
% xlabel('Date')
% ylabel('Disolved O_2')
% 
% %% Weather Data plots
% load Weather_data.mat
% 
% weather.date_time = datetime (weather.year,weather.month,weather.day,weather.hour,weather.min,0,0);
% 
% 
% figure
% plot(weather.date_time,weather.air_temp)
% xlabel('Date')
% ylabel('Air temperature (c)')
% xlim([min(Hanna.date_time) max(Hanna.date_time)])
% 
% figure
% plot(weather.date_time,weather.preasure)
% xlabel('Date')
% ylabel('Barometric Preasure (mbar)')
% xlim([min(Hanna.date_time) max(Hanna.date_time)])
% 
% figure
% plot(weather.date_time,weather.realitive_humidity)
% xlabel('Date')
% ylabel('Realitive Humidity (%)')
% xlim([min(Hanna.date_time) max(Hanna.date_time)])
% 
% figure
% plot(weather.date_time,weather.wind_speed)
% xlabel('Date')
% ylabel('Wind Speed (m/s)')
% xlim([min(Hanna.date_time) max(Hanna.date_time)])
% 

