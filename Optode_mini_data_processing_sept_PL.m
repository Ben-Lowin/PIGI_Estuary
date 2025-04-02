close all
clear all
clc

% Ben Lowin
% October 24th 2024
% calculate the NCP, PP and respiration for each incubation
% determin the YSI means for each station
% save out to a single insitu data file

% no NCP data - need to rerun from that section after processing

%% load in the incubation data
load Optode_mini_data_PL02.mat


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

day_of_sample = [ones(1,24),ones(1,24).*2,ones(1,24).*3,ones(1,24).*4, ...
    ones(1,24).*5,ones(1,24).*6,ones(1,12).*7,ones(1,24).*8,ones(1,24).*9 ...
    ,ones(1,12).*10];

% Initialize variables to store results
low_positive_means = [];
low_positive_stds = [];
low_negative_means = [];
low_negative_stds = [];
high_positive_means = [];
high_positive_stds = [];
high_negative_means = [];
high_negative_stds = [];

for day = 1:10
    %creat a day index
    day_index = find(day_of_sample == day);

    %range of indicies
    i = day_index(1):2:day_index(end-1);
    ii = day_index(2):2:day_index(end);

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
    high_negative_mean high_negative_std low_tide_i startIndex day day_index

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

% figure 
% errorbar(bottle.HT_respiration_mean,bottle.HT_respiration_std,'r*')
% hold on
% errorbar(bottle.HT_photo_mean,bottle.HT_photo_std,'g*')
% errorbar(bottle.HT_NCP_mean,bottle.HT_NCP_std,'b*')
% xlim([0 11])
% xticklabels({'', '09/16', '09/17', '09/18','09/19', '09/20', '09/21', ...
% '09/22','09/23','09/24','09/25',''});
% ylabel('mmol O_2/m^3/day')
% title('High Tide incubations')
% legend('repiration','photosynthesis','NCP')
% hold off
% 
% figure 
% errorbar(bottle.LT_respiration_mean,bottle.LT_respiration_std,'r*')
% hold on
% errorbar(bottle.LT_photo_mean,bottle.LT_photo_std,'g*')
% errorbar(bottle.LT_NCP_mean,bottle.LT_NCP_std,'b*')
% xlim([0 11])
% xticklabels({'', '09/16', '09/17', '09/18','09/19', '09/20', '09/21', ...
% '09/22','09/23','09/24','09/25',''});
% ylabel('mmol O_2/m^3/day')
% title('Low Tide incubations')
% legend('repiration','photosynthesis','NCP')
% hold off


%% put date time in bottles so they can be plotted

bottle.sample_dates(1) = datetime(2024,9,16,7,29,00,00);
bottle.sample_dates(2) = datetime(2024,9,16,13,49,00,00);
bottle.sample_dates(3) = datetime(2024,9,17,8,22,00,00);
bottle.sample_dates(4) = datetime(2024,9,17,14,58,00,00);
bottle.sample_dates(5) = datetime(2024,9,18,9,14,00,00);
bottle.sample_dates(6) = datetime(2024,9,18,15,52,00,00);
bottle.sample_dates(7) = datetime(2024,9,19,10,06,00,00);
bottle.sample_dates(8) = datetime(2024,9,19,16,46,00,00);
bottle.sample_dates(9) = datetime(2024,9,20,10,58,00,00);
bottle.sample_dates(10) = datetime(2024,9,20,17,38,00,00);
bottle.sample_dates(11) = datetime(2024,9,21,11,52,00,00);
bottle.sample_dates(12) = datetime(2024,9,21,18,30,00,00);
bottle.sample_dates(13) = datetime(2024,9,22,12,49,00,00);
bottle.sample_dates(14) = datetime(2024,9,23,7,35,00,00);
bottle.sample_dates(15) = datetime(2024,9,23,13,49,00,00);
bottle.sample_dates(16) = datetime(2024,9,24,8,33,00,00);
bottle.sample_dates(17) = datetime(2024,9,24,14,50,00,00);
bottle.sample_dates(18) = datetime(2024,9,25,9,35,00,00);

% put NCP into a single varable based on the above datetime
count = 0;
for i = 1:length(bottle.HT_NCP_mean)

      if isnan (bottle.HT_NCP_mean(i))

      else
        count = count + 1;
        bottle.all_NCP(count) = bottle.HT_NCP_mean(i);
        bottle.all_NCP_std(count)  = bottle.HT_NCP_std(i);
        bottle.all_photo(count) = bottle.HT_photo_mean(i);
        bottle.all_photo_std(count) = bottle.HT_photo_std(i);
        bottle.all_resp(count) = bottle.HT_respiration_mean(i);
        bottle.all_resp_std(count) = bottle.HT_respiration_std(i);

     end
    
    if isnan (bottle.LT_NCP_mean(i))
    
    else
        count = count + 1;
        bottle.all_NCP(count) = bottle.LT_NCP_mean(i);
        bottle.all_NCP_std(count)  = bottle.LT_NCP_std(i);
        bottle.all_photo(count) = bottle.LT_photo_mean(i);
        bottle.all_photo_std(count) = bottle.LT_photo_std(i);
        bottle.all_resp(count) = bottle.LT_respiration_mean(i);
        bottle.all_resp_std(count) = bottle.LT_respiration_std(i);
    end
    
end

clear count i 

%% load NCP for comparison


load PL02_NCP_01.mat

% figure
% plot(pigi_dat.datetime,pigi_dat.NCP)
% xlabel('Date')
% ylabel('NCP (umol O_2/m^2/day)')
% 
% figure
% plot(pigi_dat.datetime,pigi_dat.n2_calc.do2n2)
% xlabel('Date')
% ylabel('delta O_2/N_2 (uM O_2/day)')


%%

NCP_m3=pigi_dat.NCP./3.5;



for i = 1:18

range = find (pigi_dat.datetime < bottle.sample_dates(i) + hours(1) & ...
    pigi_dat.datetime > bottle.sample_dates(i) - hours(1));

NCP_avr(i)=mean(NCP_m3(range),"omitnan");
NCP_std(i)=std(NCP_m3(range),"omitnan");

clear range
end


% figure
% plot(pigi_dat.datetime,NCP_m3)
% xlabel('Date')
% ylabel('NCP (mmol O_2/m^3/day)')

%%

% figure
% plot(1:18, NCP_avr,'*')
% errorbar(1:18,NCP_avr,NCP_std,NCP_std,'*')
% hold on
% errorbar(1:18,bottle.all_NCP,bottle.all_NCP_std,bottle.all_NCP_std,'*')
% xlim([0 19])
% ylabel('NCP pigi (mM O_2/m^3/day)')
% legend ('Pigi data','Bottle data')


difference_bottle_pigi = NCP_avr - bottle.all_NCP;

bottle.pigi_ncp_avr = NCP_avr;
bottle.pigi_ncp_std = NCP_std;
bottle.pigi_NCP_m3 = NCP_m3;
bottle.pigi_datetime = pigi_dat.datetime;

% save('NCP_bottle_data_PL02','bottle')
