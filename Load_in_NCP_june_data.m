close all
clear all
clc

% Ben Lowin
% June 26 2024
% load all the data into matlab

%% Load in the weather data

weather_inital = readtable("weather_may_june_data_long.csv");
%date and time
weather.year = weather_inital.YearUnits;
weather.month = weather_inital.MonthUnits;
weather.day = weather_inital.DayUnits;
weather.hour = weather_inital.HourUnits;
weather.min = weather_inital.MinuteUnits;
%weather data
weather.wind_speed = weather_inital.RelativeWindSpeedM_s;
weather.air_temp = weather_inital.TemperatureC;
weather.realitive_humidity = weather_inital.RelativeHumidity_RH;
weather.preasure = weather_inital.PressureHectoPascals;
weather.units = ["wind speed - m/s","air temp - c", "realitive  humidity - %",...
    "preasure - mbar"];

clear weather_inital

%% Load in the YSI data

files = dir('C:\Users\BenL\Documents\MATLAB\PIGI\June_dock_NCP\YSI_data');

% List of data names
data_names = ["monday_high_in", "monday_low_in", "monday_high_out", "tuesday_high_in",...
    "monday_low_out", "tuesday_low_in", "tuesday_high_out", "wednesday_high_in", ...
    "tuesday_low_out", "wednesday_low_in", "wednesday_high_out", "thursday_high_in", ...
    "wednesday_low_out", "thursday_low_in", "thursday_high_out", "friday_high_in", ...
    "thursday_low_out", "friday_low_in"];

% Initialize structure
YSI = struct();

% Loop through each file and assign data to corresponding fields
for idx = 1:18
    holder = readtable(fullfile(files(idx+2).folder, files(idx+2).name)); % Read the table
    
    data_name = data_names(idx); % Current data name
    YSI.(data_name).date = holder.Date;
    YSI.(data_name).time = holder.Time;
    YSI.(data_name).temp = holder.x_C;
    YSI.(data_name).disolved_oxygen = holder.DO_RTB;
    YSI.(data_name).salinity = holder.SAL_PSU;
    YSI.(data_name).pH = holder.pH;
    YSI.(data_name).chl = holder.ChlUg_L;
end
YSI.units = ["temperature - c", "disolved oxygen - %", "salinity - psu", ...
    "pH - pH", "chlorophyll - ug/l"];

clear files data_name data_names idx

%% Load in Hanna Data

hanna_initial = readtable("Hanna_june_data.xls", "Sheet"," Log data - 1");

Hanna.date = hanna_initial.Date;
Hanna.time = hanna_initial.Time;
Hanna.temp = hanna_initial.Temp___C_;
Hanna.pH = hanna_initial.pH;
Hanna.salinity = hanna_initial.Sal__psu_;
Hanna.density = hanna_initial.SigmaT_sT_+1000;
Hanna.Disolved_oxygen = hanna_initial.D_O__ppm_;
Hanna.units =["temperature - c", "pH - pH", "salinity - psu", ...
    "density - g/kg", "disolved oxygen - ppm"];

clear hanna_initial

%% Load Mini Optode data
% Get list of files
files = dir('C:\Users\BenL\Documents\MATLAB\PIGI\June_dock_NCP\MINI_optode_data\Master\*.txt');

% Initialize a structure to store the data
Optode_mini = struct();

% Loop through each file
for idx = 1:length(files)
    file = fullfile(files(idx).folder, files(idx).name); % Get the full path of the file
    opts = detectImportOptions(file, 'Delimiter', ';'); % Create import options with the delimiter set to semicolon
    opts.VariableNamesLine = 58;
    opts.DataLines = 59;
    holder = readtable(file, opts); % Read the table using the specified options

    % Extract the date and sample information from the filename
    name = files(idx).name;
    dateField = ['d' name(1:8)];
    sampleField = name(10:end-4); % Remove the file extension
    
    % Store the data in the structure
    Optode_mini.(dateField).(sampleField).oxygen = holder.oxygen_umol_L;
    Optode_mini.(dateField).(sampleField).temperature = holder.temp__C;
    Optode_mini.(dateField).(sampleField).amp = holder.amp;
    
end
Optode_mini.units = ["oxygen - umol/l", "temperature - c", "amp - unitless"];

clear dateField file files holder idx name opts sampleField
clc

%% Load water level data

water_initial = readtable("water_level_PL_june.csv");

% Convert epoch time from milliseconds to seconds
water_epochTimeSeconds = water_initial.Timestamp / 1000;

% Convert epoch time to datetime
water_dateTime = datetime(water_epochTimeSeconds, 'ConvertFrom', 'posixtime');

water_height.time = water_dateTime;
water_height.depth = water_initial.DistanceMeters;

clear water_initial water_epochTimeSeconds water_dateTime

%% Save out the varables

% cd C:\Users\BenL\Documents\MATLAB\PIGI\June_dock_NCP
% 
% save('Hanna_data.mat','Hanna')
% save('Optode_mini_data.mat',"Optode_mini")
% save('Weather_data.mat',"weather")
% save("YSI_data.mat","YSI")
% save('Water_height',"water_height")
