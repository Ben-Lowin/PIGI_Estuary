close all
clear all
clc

% Playing with June_dock data and Net ecosystem production
% based on methods from Wang et al 2017
% Inorganic carbon and oxygen dynamics in a marsh-dominated estuary

%looking to deploy the day/night oxygen change assumption

%% load in data
load PL01_NCP_03.mat

load NCP_bottle_data.mat

sun_time = readtable("sun_rise_set_june_2024.xlsx");

 %% find the o2 for sunrise and sun set

sunrise = datetime(sun_time.sunRise);
sunset = datetime(sun_time.sunSet);

range = 360; %is one hour and the same range as used in the other data processing

for in = 1:length(sunrise)

    [M,I]=min(abs(pigi_dat.datetime-sunrise(in)));
    if M<hours(1)
        sunrise_m(in)= mean(pigi_dat.ts_cor.o2uM(I-range:I+range),'omitnan');
        sunrise_std(in)= std(pigi_dat.ts_cor.o2uM(I-range:I+range),'omitnan');
        sunrise_t(in) = sunrise(in);
        sunrise_flux(in) = mean(pigi_dat.flux_02(I-range:I+range),'omitnan');
        sunrise_flux_std(in) = std(pigi_dat.flux_02(I-range:I+range),'omitnan');

    end

    [MM,II]=min(abs(pigi_dat.datetime-sunset(in)));
    if MM<hours(1)
        sunset_m(in)= mean(pigi_dat.ts_cor.o2uM(II-range:II+range),'omitnan');
        sunset_std(in)= std(pigi_dat.ts_cor.o2uM(II-range:II+range),'omitnan');
        sunset_t(in) = sunset(in);
        sunset_flux(in) = mean(pigi_dat.flux_02(II-range:II+range),'omitnan');
        sunset_flux_std(in) = std(pigi_dat.flux_02(II-range:II+range),'omitnan');
    end

end

clear I II in M MM range suntime

%% plot o2 to check it matches

figure

% Left Y-axis: Oxygen, sunrise, sunset
yyaxis left
plot(pigi_dat.datetime, pigi_dat.ts_cor.o2uM, 'b-')
ylabel('Oxygen (\mumol/L or mmol/m^3)')
hold on
plot(sunrise_t, sunrise_m, 'ro')
plot(sunset_t, sunset_m, 'ko')

% Right Y-axis: Tide depth
yyaxis right
plot(pigi_dat.datetime, pigi_dat.tide_depth, 'k-')
ylabel('Tide Depth')
ylim([0 12])

% Customize axes, legend, etc.
legend('Oxygen','Sunrise','Sunset','Tide Depth','Location','best')
xlabel('Date/Time')
hold off

figure 
plot(pigi_dat.datetime,pigi_dat.flux_02)
hold on
plot(sunrise_t, sunrise_flux, 'ro')
plot(sunset_t, sunset_flux, 'ko')
ylabel('diffusive flux (mmol/m^2/hr)')
xlabel('Date/Time')
hold off

%%

%how many hours in the respiration calculation

respiration_hours = sunrise_t(3:end) - sunset_t(2:end);
respiration_hours_simple = hours(respiration_hours);
flux_respiration = (sunrise_flux(3:end) + sunset_flux(2:end))./2.*respiration_hours_simple;
flux_respiration_std = sqrt(sunrise_flux_std(3:end).^2 + sunset_flux_std(2:end).^2)./2.*respiration_hours_simple;

respiration = (sunrise_m(3:end)-sunset_m(2:end)- flux_respiration)./mean(respiration_hours_simple);
respiration_sdt = sqrt (sunrise_std(3:end).^2+sunset_std(2:end).^2+flux_respiration_std.^2)./mean(respiration_hours_simple);

NDP_hours = sunset_t(2:end) - sunrise_t(2:end-1);
NDP_hours_simple = hours(NDP_hours);
flux_NDP= (sunset_flux(2:end) + sunrise_flux(2:end-1))./2.*NDP_hours_simple;
flux_NDP_std = sqrt(sunset_flux_std(2:end).^2 + sunrise_flux_std(2:end-1).^2)./2.*NDP_hours_simple;

NDP =(sunset_m(2:end)- sunrise_m(2:end-1) - flux_NDP)./mean(NDP_hours_simple);
NDP_std = sqrt (sunset_std(2:end).^2 + sunrise_std(2:end-1).^2 + flux_NDP_std.^2)./mean(NDP_hours_simple);


GPP = NDP-respiration;
GPP_std = sqrt (NDP_std.^2 + respiration_sdt.^2);

%I need a depth value to divide out
depth = mean(pigi_dat.tide_depth,'omitnan');

NEP = (GPP.*mean(NDP_hours_simple) + respiration.*24) ./depth  ;
NEP_std = sqrt((mean(NDP_hours_simple) .* GPP_std).^2 + (24 .* respiration_sdt).^2)./depth;




%%
figure
errorbar(bottle.sample_dates,bottle.all_NCP,bottle.all_NCP_std)
hold on
errorbar(bottle.sample_dates,bottle.pigi_ncp_avr.*3.5./depth,bottle.pigi_ncp_std*3.5./depth)
errorbar(sunrise_t(2:6),NEP,NEP_std)

legend('Bottle NCP', 'PIGI NCP','Odum NCP')
xlabel('Date')
ylabel('NCP (mM/m^2/day)')
hold off

%%

PL_01_three.btl_time = bottle.sample_dates;
PL_01_three.btl_ncp = bottle.all_NCP;
PL_01_three.btl_std = bottle.all_NCP_std;
PL_01_three.pigi_time = bottle.sample_dates;
PL_01_three.pigi_ncp = bottle.pigi_ncp_avr.*3.5./depth;
PL_01_three.pigi_std = bottle.pigi_ncp_std*3.5./depth;
PL_01_three.odum_time = sunrise_t(2:6);
PL_01_three.odum_ncp = NEP;
PL_01_three.odum_std = NEP_std;


%save('PL01_three_plot', 'PL_01_three')




%%
% figure
% plot(pigi_dat.datetime,pigi_dat.ts_cor.windspeed)
% xlabel('time')
% ylabel('Wind speed (m/s)')
% xlim([pigi_dat.datetime(5124) pigi_dat.datetime(48257)])
% 
% %% sum of sunlight per day
% 
% for in = 1:length(sunrise_t)-1
% 
%     aa = find (pigi_dat.datetime > sunrise_t(in) & pigi_dat.datetime < sunrise_t(in+1));
% 
%     solar_mean(in) = mean(pigi_dat.ts_cor.solar_rad(aa),'omitnan');
% 
% end
% 
% figure
% plot(sunrise_t(1:end-1),solar_mean)
% xlabel('time')
% ylabel('mean solar radiation (W/M^2)')
% 
% figure 
% plot(pigi_dat.datetime, pigi_dat.ts_cor.sun_hours)
% xlabel('time')
% ylabel('solar radiation (W/M^2)')
% xlim([pigi_dat.datetime(5124) pigi_dat.datetime(48257)])