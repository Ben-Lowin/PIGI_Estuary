close all
clear all
clc

% Ben Lowin
% October 24th 2024
% load all the data into matlab

%% Load in the weather data

weather_inital = readtable("Historic_weather_data_PL_sept_aug.csv");
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


%% Load in Hanna Data

hanna_initial = readtable("Hanna_sept_data.xlsx");

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
files = dir('C:\Users\BenL\Documents\MATLAB\PIGI\Sep_dock_NCP\Mini_optode_data\Master\*.txt');

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

water_initial = readtable("Waterlevel_PL_sept.csv");

% Convert epoch time from milliseconds to seconds
water_epochTimeSeconds = water_initial.Timestamp / 1000;

% Convert epoch time to datetime
water_dateTime = datetime(water_epochTimeSeconds, 'ConvertFrom', 'posixtime');

water_height.time = water_dateTime;
water_height.depth = water_initial.DistanceMeters;

clear water_initial water_epochTimeSeconds water_dateTime

%% Save out the varables

% cd C:\Users\BenL\Documents\MATLAB\PIGI\Sep_dock_NCP
% 
% save('Hanna_data_PL02.mat','Hanna')
% save('Optode_mini_data_PL02.mat',"Optode_mini")
% save('Weather_data_PL02.mat',"weather")
% save('Water_height_PL02',"water_height")
