close all
clear all
clc

% Playing with sep_dock data and Net ecosystem production
% based on methods from Wang et al 2017
% Inorganic carbon and oxygen dynamics in a marsh-dominated estuary

%looking to deploy the day/night oxygen change assumption

%%
load PL03_NCP_03.mat
load NCP_bottle_data_PL03.mat
sun_time = readtable("december_sunrise_PL.xlsx");

 %% find the o2 for sunrise and sun set

sunrise = datetime(sun_time.sunrise);
sunset = datetime(sun_time.sunset);

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
plot(pigi_dat.datetime,pigi_dat.ts_cor.o2uM)
ylabel('Oxygen (umol/L)')
hold on 

plot(sunrise_t,sunrise_m,'ro')
plot(sunset_t,sunset_m,'ko')
legend('Oxygen','sun rise','sun set')
hold off

%%
%how many hours in the respiration calculation

respiration_hours = sunrise_t(3:end) - sunset_t(2:end-1);
respiration_hours_simple = hours(respiration_hours);
flux_respiration = (sunrise_flux(3:end) + sunset_flux(2:end-1))./2.*respiration_hours_simple;
flux_respiration_std = sqrt(sunrise_flux_std(3:end).^2 + sunset_flux_std(2:end-1).^2)./2.*respiration_hours_simple;

respiration = (sunrise_m(3:end)-sunset_m(2:end-1)- flux_respiration)./mean(respiration_hours_simple);
respiration_sdt = sqrt (sunrise_std(3:end).^2+sunset_std(2:end-1).^2+flux_respiration_std.^2)./mean(respiration_hours_simple);

NDP_hours = sunset_t(2:end-1) - sunrise_t(2:end-1);
NDP_hours_simple = hours(NDP_hours);
flux_NDP= (sunset_flux(2:end-1) + sunrise_flux(2:end-1))./2.*NDP_hours_simple;
flux_NDP_std = sqrt(sunset_flux_std(2:end-1).^2 + sunrise_flux_std(2:end-1).^2)./2.*NDP_hours_simple;

NDP =(sunset_m(2:end-1)- sunrise_m(2:end-1) - flux_NDP)./mean(NDP_hours_simple);
NDP_std = sqrt (sunset_std(2:end-1).^2 + sunrise_std(2:end-1).^2 + flux_NDP_std.^2)./mean(NDP_hours_simple);


GPP = NDP-respiration;
GPP_std = sqrt (NDP_std.^2 + respiration_sdt.^2);

%I need a depth value to divide out
depth = mean(pigi_dat.tide_depth,'omitnan');

NEP = (GPP.*mean(NDP_hours_simple) + respiration.*24) ./depth  ;
NEP_std = sqrt((mean(NDP_hours_simple) .* GPP_std).^2 + (24 .* respiration_sdt).^2)./depth;



%% plotting all three


figure
errorbar(bottle.sample_dates,bottle.all_NCP,bottle.all_NCP_std)
hold on
errorbar(bottle.sample_dates,bottle.pigi_ncp_avr.*3.5./depth,bottle.pigi_ncp_std*3.5./depth)
errorbar(sunrise_t(2:5),NEP,NEP_std)

legend('Bottle NCP', 'PIGI NCP','Odum NCP')
xlabel('Date')
ylabel('NCP (mM/m^2/day)')
hold off

%%

PL_03_three.btl_time = bottle.sample_dates;
PL_03_three.btl_ncp = bottle.all_NCP;
PL_03_three.btl_std = bottle.all_NCP_std;
PL_03_three.pigi_time = bottle.sample_dates;
PL_03_three.pigi_ncp = bottle.pigi_ncp_avr.*3.5./depth;
PL_03_three.pigi_std = bottle.pigi_ncp_std*3.5./depth;
PL_03_three.odum_time = sunrise_t(2:5);
PL_03_three.odum_ncp = NEP;
PL_03_three.odum_std = NEP_std;


%save('PL03_three_plot', 'PL_03_three')


% %% potential explaitory varables?
% 
% figure
% plot(pigi_dat.datetime,pigi_dat.ts_cor.windspeed)
% xlabel('date')
% ylabel('Wind speed (m/s)')
% % windspeed seems to plot well to justify the change
% % pH has not trend
% % Air temp does aligh but is not as strong as wind
% % solar radation does not match too well, might need to sum per day
% 
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
% 
% figure 
% plot(pigi_dat.datetime, pigi_dat.ts_cor.solar_rad)
% xlabel('time')
% ylabel('solar radiation (W/M^2)')
% xlim ([pigi_dat.datetime(5370) pigi_dat.datetime(39884) ])
% 
% figure 
% plot(pigi_dat.datetime, pigi_dat.ts_cor.windspeed)
% xlabel('time')
% ylabel('Wind speed (m/s)')
% xlim ([pigi_dat.datetime(5370) pigi_dat.datetime(39884) ])
% 
% %%

