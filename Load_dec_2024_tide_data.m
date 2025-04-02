

%load in the december2024_tide_PL data
%note load by hand
%%
inp = december2024tidePL;

dateStrings = inp.Date;
timeStrings = inp.Day;

% Convert categorical dates to cell array and then to datetime
dateCells = cellstr(dateStrings);
dates = datetime(dateCells, 'InputFormat', 'yyyy/MM/dd');

% Convert time strings to datetime (time-only)
times = datetime(timeStrings, 'InputFormat', 'hh:mm a');

% Combine date and time
dateTimeCombined = dates + timeofday(times);
depth = inp.Time;

%%

raw_date_time = dateTimeCombined(1) : seconds(10) : dateTimeCombined(end);


interpolated_depth=interp1(dateTimeCombined,depth,raw_date_time,'spline');

figure
plot(raw_date_time,interpolated_depth)

tide.time = raw_date_time;
tide.depth = interpolated_depth;

save('NOAA_dec_tide','tide')
