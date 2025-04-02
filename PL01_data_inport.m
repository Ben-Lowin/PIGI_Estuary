close all
clear all
clc

%Ben Lowin
%06/26/2024
%import data from Load_in_NCP_june_data into a PIGI readable

%% underway data file directory
cd ('C:\Users\BenL\Documents\MATLAB\PIGI\June_dock_NCP')

%% Read in cruise data_ full
% translate to a CSV using excel
%change for new data
load("Hanna_data.mat")
load("Weather_data_PL01_02.mat")

%% make a date time vector for each

date_h = Hanna.date + Hanna.time;
date_w = datetime(weather.year, weather.month, weather.day, weather.hour, ...
    weather.min,zeros(length(weather.day),1),zeros(length(weather.day),1) );

index_weather_match = find (date_w >= datetime(2024, 06, 16, 9, 40, 0, 0));


%% interpolate the data to be the same time step (1/min)
% figures
% plot(date_h,Hanna.salinity) 
% 
% figure 
% plot(date_h, Hanna.Disolved_oxygen) 
% 
% figure
% plot(date_h, Hanna.temp) %- I am not sure I trust the declining water temp
%turns out that the water temp is good. there is about 0.2 c disagrement
%between the YSI and Hanna - this might be depth of deploymnet? or
%calibration?

interp_grid = datetime(2024, 06, 16, 9, 40, 0, 0) :minutes(1): datetime(2024, 06, 22, 15, 20, 0, 0);
%grid on which i will pull data points

a =interp1(date_h, Hanna.temp,interp_grid);

% might want to apply smothing (moving avrage or box car)
% figure
% plot(interp_grid,a)
% hold on
% plot(date_h, Hanna.temp)
%%
Inported_data = struct();

Inported_data.time=interp_grid;
Inported_data.temperature=interp1(date_h,Hanna.temp,interp_grid);
Inported_data.salinity=interp1(date_h,Hanna.salinity,interp_grid);
Inported_data.density=interp1(date_h,Hanna.density,interp_grid);
Inported_data.disolved_oxygen=interp1(date_h,Hanna.Disolved_oxygen,interp_grid);
Inported_data.pH =interp1(date_h,Hanna.pH,interp_grid);

Inported_data.wind_speed=interp1(date_w,weather.wind_speed,interp_grid);
Inported_data.barometeric =interp1(date_w,weather.preasure,interp_grid);
Inported_data.air_temp=interp1(date_w,weather.air_temp,interp_grid);
Inported_data.realitive_humidity=interp1(date_w,weather.realitive_humidity,interp_grid);
Inported_data.sun_hours = interp1(date_w, weather.solar_radation, interp_grid);
Inported_data.solar_rad = interp1(date_w, weather.solar_rad, interp_grid);

% do I need a true depth at location?

Inported_data.units=["time - datenum", "temperature - c", "salinity - PSU", ...
    "desnity - g/kg", "Disolved Oxygen - ppm", "pH - pH","wind speed - m/s", ...
    "barometic preasure - mbar", "air temperature - c", "realitive humidity - %",...
    "Solar hours - hours","Solar radation - W/m^2"];
    

%% save out
%change for new data

% save('PL01_master',"Inported_data")