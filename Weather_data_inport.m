weather.air_temp = weathermayjunedatalong.TemperatureC;
weather.preasure = weathermayjunedatalong.PressureHectoPascals;
weather.realitive_humidity = weathermayjunedatalong.RelativeHumidityRH;
weather.wind_speed = weathermayjunedatalong.RelativeWindSpeedMs;
weather.year = weathermayjunedatalong.YearUnits;
weather.month = weathermayjunedatalong.MonthUnits;
weather.day = weathermayjunedatalong.DayUnits;
weather.hour = weathermayjunedatalong.HourUnits;
weather.min = weathermayjunedatalong.MinuteUnits;
weather.solar_radation = weathermayjunedatalong.SolarRadiationWm2;
weather.solar_rad = weathermayjunedatalong.SunshineHoursHours;

weather.datetime = datetime(weather.year, weather.month, weather.day, weather.hour, weather.min, 0);

save('weather_data_PL01_02',"weather")