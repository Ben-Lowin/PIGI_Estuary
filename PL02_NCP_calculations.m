close all
clear all
clc

%% Load data
%Load the processed pigi data
load PL_02_deployment_PIGI_RAW_v20250212.mat
%load the imported temp/salintity data
load("PL02_master.mat")
%load historical wind speed data
load("Weather_data_PL02.mat")
%load water heigh
load("Water_height_PL02.mat")


%% set the temp, salinity and wind to be the right points
%change the output to be in the right units
%and change the orientation
do2n2=pigi_dat.n2_calc.do2n2'./100;
S=pigi_dat.n2_calc.sal';
T=pigi_dat.n2_calc.sst';

% interpolate the wind speed to the PIGI data
u10mat_obs=interp1(Inported_data.time, Inported_data.wind_speed, pigi_dat.datetime)';

%% calculate the water depth (mld)
% assuming that the zero level is mean tide level
%assuming that this is 1.5m above lower low tide (chart values)
%assuming that PL doc is at 1m depth at lower low tide (per chart value)

PL_water_depth = water_height.depth +1.9;
% plot(water_height.time,PL_water_depth)
% Looks good

% Ensure unique and sorted values for water_height.time
[water_time_unique, idx] = unique(water_height.time);

% Use the indices to get corresponding PL_water_depth values
PL_water_depth_unique = PL_water_depth(idx);

tide_time = water_time_unique - hours(3.916);

% Now perform the interpolation with unique sample points
mld = interp1(tide_time, PL_water_depth_unique, pigi_dat.datetime)';

clear PL_water_depth
%% Make historical wind speed 

%define varable
ws = weather.wind_speed;

date_w = datetime(weather.year, weather.month, weather.day, weather.hour, ...
    weather.min,zeros(length(weather.day),1),zeros(length(weather.day),1) );

% Ensure unique and sorted values for date_w
[date_w_unique, idx_w] = unique(date_w);

% Use the indices to get corresponding ws (wind speed) values
ws_unique = ws(idx_w);

% Now perform the interpolation with unique values
pigi_wind = interp1(date_w_unique, ws_unique, pigi_dat.datetime)';


%define matrix lenghts
N=length(pigi_wind);
M=121;

% Create a timetable
windTimetable = timetable(date_w, ws);

% Average data to 6-hour intervals
averagedWindTimetable = retime(windTimetable, 'regular', 'mean', 'TimeStep', hours(6));

u10mat_new = NaN(N,M);

for idx = 1:N
    if ~isnat(pigi_dat.datetime(idx))
        %finds the point to start at (based on the current data
        [i, d] = min (abs(pigi_dat.datetime(idx) -averagedWindTimetable.date_w));
        %sets the old data into place (for the last 23 6 hour avrages)
        u10mat_new(idx,1:M-1) = averagedWindTimetable.ws(d-120:d-1);
        % adds the new data into the last row
        u10mat_new(idx,M) =pigi_wind(idx);
    else
        u10mat_new(idx,:) = nan;
    end
end

clear idx averagedWindTimetable N M pigi_wind ws date_w i d windTimetable

%% Define last input varables

% INPUT:
%       do2n2,S,T,mld
%                 = arrays for delO2/N2 [decimal %], sea surface salinity
%                     [PSU], sea surface temperature [C] and mixed layer
%                     depth [m]. Note that all arrays should be the same 
%                     size [1xN].
%       u10mat    = wind speed (at 10m) matrix [m/s] of size [MxN] where M
%                     corresponds with the number of observations backwards 
%                     in time. M = Ndays / dt + 1;
%                     This matrix is used to calculate a weighted O2 piston velocity.
%             NOTE: u10mat(:,M) = wind speed corresponding with time of
%             underway observations. u10mat(:,1) = most "historic" wind
%             speed.
%       Ndays     = number of days over which to perform piston velocity
%                     weighting (recommend 30 or 60)
%       dt        = time increment, in days, of u10mat and mldmat
%Days before work of data
Ndays=30;
%time step
dt=.25;
%run mode
mod="w14";
%M=Ndays/(dt+1) = 24


%% run the NCP calculation
% [ncp,kw_o2,tro2,wt_t] = calc_o2n2_ncp(do2n2,S,T,mld,u10mat_new,Ndays,dt,mod);

[ncp,kw_o2,tro2,wt_t] = calc_o2n2_ncp_estuary(do2n2,S,T,mld,u10mat_new,Ndays,dt,mod);

%% Plot O2 and N2 vs NCP 
%--- Plot
tt= pigi_dat.datetime;
% 
% figure
%     subplot(8,1,1); hold on
%         plot(tt,do2n2*100,'k','linewidth',2)
%         ylabel('\DeltaO2/N2 [%]')
%     subplot(8,1,2); hold on
%         plot(tt,ncp,'r','linewidth',2)
%         ylabel('NCP [mmol O2/m2/d]')
%         set(gca,'yaxisloc','right')
%     subplot(8,1,3); hold on
%         plot(tt,T,'k','linewidth',2)
%         ylabel('SST')
%     subplot(8,1,4); hold on
%         plot(tt,S,'k','linewidth',2)
%         ylabel('Sal.')
%         set(gca,'yaxisloc','right')
%     subplot(8,1,5); hold on
%         plot(tt,u10mat_new(:,end),'k','linewidth',2)
%         plot(tt,u10mat_new);
%         plot(tt,u10mat_new(:,end),'k','linewidth',2)
%         ylabel('u10 [m/s]')
% %         legend('at time of underway obs.','during weighting period','location','eastoutside')
%     subplot(8,1,6); hold on
%         plot(tt,mld,'k','linewidth',2)
%         ylabel('MLD [m]')
%         set(gca,'yaxisloc','right')
%         axis ij
%     subplot(8,1,7); hold on
%         plot(tt,tro2,'k','linewidth',2)
%         plot(tt,wt_t,'b','linewidth',2)
%         ylabel('Time [days]')
%         legend('O2 residence time','Adjusted weighting time','location','eastoutside')
%     subplot(8,1,8); hold on
%         plot(tt,kw_o2,'k','linewidth',2)
%         ylabel('ko2 [m/d]')
%         set(gca,'yaxisloc','right')
%         xlabel('Time')



%%
figure
subplot(3,1,1)
        plot(tt,do2n2*100,'k','linewidth',2)
        ylabel('\DeltaO2/N2 [%]')
subplot(3,1,2)
        plot(tt,ncp,'r','linewidth',2)
        ylabel('NCP [mmol O2/m2/d]')
subplot(3,1,3)
        plot(tt,kw_o2,'k','linewidth',2)
        ylabel('ko2 [m/d]')

%% estimating oxygen diffusion flux

C_measured = pigi_dat.qc_cal.o2uM;
sat = pigi_dat.qc_cal.o2sat./100;

% 1- determin O2 at 100%
C_sat = C_measured./sat;

% 2- calucalte delta oxygen
C_delta = C_sat - C_measured;

% 3- convert K to /h
Kwh = kw_o2./24;

% 4- convert delta oxygen to (mol/m^3)
C_delta_m_cubed = C_delta.*10^-3;

% 5- Calculate flux (mmmol/m^2/hr)
Flux_o2 = Kwh.*C_delta_m_cubed*1000;

figure
plot(Flux_o2,'.')



%% adding the NCP data to the data structure and saving out

pigi_dat.tide_depth = mld;
pigi_dat.flux_02_units = 'mmol/m^2/hr';
pigi_dat.flux_02 = Flux_o2;
pigi_dat.NCP_units = 'mmol/m^2/day';
pigi_dat.NCP=ncp;
% save("PL02_NCP_03", "pigi_dat")





