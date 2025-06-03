%% load past TWL forecasts
% 
MySQLDBConnect
%tampa bay region
region = 'TBW';
%roughly position of camera
% latitude  = 27.7958;
% longitude =-82.7956;

% cam.lat = 27.796216949206798;
% cam.lon = -82.796102542950635;
latitude = 27.796216949206798;
longitude = -82.796102542950635;

%date
dateTime='2017-06-03 08:00:00';

% this will load all forecasts within the ~104 hr window before/after dateTime
[twl_output] = MySQLDBGetTWLTimeseries(longitude, latitude, dateTime, region);

% this will load all forecasts for dateTime ONLY (12 timepoints?)
% [twl_output] = MySQLDBGetTWLForecast(longitude, latitude, dateTime, region)


for ii = 1:length(twl_output)
year(ii) = str2num(twl_output(ii).forecastTime(1:4));
month(ii) = str2num(twl_output(ii).forecastTime(6:7));
day(ii) = str2num(twl_output(ii).forecastTime(9:10));
hr(ii) = str2num(twl_output(ii).forecastTime(12:13));
MM = 0;
SS = 0;
date(ii) = datenum(year(ii),month(ii),day(ii),hr(ii),MM,SS);
twl_output(ii).tide = twl_output(ii).twl-twl_output(ii).runup;
end
% only use this if wanting to compare "forecastBegins" times:
for ii = 1:length(twl_output)
year(ii) = str2num(twl_output(ii).forecastBegins(1:4));
month(ii) = str2num(twl_output(ii).forecastBegins(6:7));
day(ii) = str2num(twl_output(ii).forecastBegins(9:10));
hr(ii) = str2num(twl_output(ii).forecastBegins(12:13));
MM = 0;
SS = 0;
date(ii) = datenum(year(ii),month(ii),day(ii),hr(ii),MM,SS);
twl_output(ii).tide = twl_output(ii).twl-twl_output(ii).runup;
end

figure
plot(date,vertcat(twl_output.twl))
hold on
plot(date,vertcat(twl_output.runup))
plot(date,vertcat(twl_output.twl)-vertcat(twl_output.runup))
datetick('x',6)
xlabel('date');ylabel('elevation')
legend('TWL','Runup','Tide')
title('TWL Forecast')
% h=fill([date,fliplr(date)],[vertcat(twl_output.twl95);flipud(vertcat(twl_output.twl05))],...
%     [0 0.8 0.9],'edgecolor',[0 0.8 0.9]);
% uistack(h,'bottom')
% h=fill([date,fliplr(date)],[vertcat(twl_output.runup95);flipud(vertcat(twl_output.runup05))],...
%     [1 0.6 0.6],'edgecolor',[1 0.6 0.6]);
% uistack(h,'bottom')


%% download any stretch of data by day with 'most recent' forecast
day1=datenum(2017,5,2);
day2=datenum(2017,5,5);
clear twl_output1
for jj=day1:day2
    dateTime = datestr(jj,'yyyy-mm-dd HH:MM:SS')
    [twl_output] = MySQLDBGetTWLTimeseries(longitude, latitude, dateTime, region);
    
    
    for ii = 1:length(twl_output)
        year(ii) = str2num(twl_output(ii).forecastTime(1:4));
        month(ii) = str2num(twl_output(ii).forecastTime(6:7));
        day(ii) = str2num(twl_output(ii).forecastTime(9:10));
        hr(ii) = str2num(twl_output(ii).forecastTime(12:13));
        MM = 0;
        SS = 0;
        date(ii) = datenum(year(ii),month(ii),day(ii),hr(ii),MM,SS);
        twl_output(ii).tide = twl_output(ii).twl-twl_output(ii).runup;
    end
    
    a=find(date(1,:)==jj);
    b=find(date(1,:)==jj+1); b=b-1;
    if jj==day1
        twl_output1=twl_output(a:b);
    end
    twl_output1=[twl_output1;twl_output(a:b)];
        
    clear year month day hr date
end

twl_output=twl_output1;
    for ii = 1:length(twl_output)
        year(ii) = str2num(twl_output(ii).forecastTime(1:4));
        month(ii) = str2num(twl_output(ii).forecastTime(6:7));
        day(ii) = str2num(twl_output(ii).forecastTime(9:10));
        hr(ii) = str2num(twl_output(ii).forecastTime(12:13));
        MM = 0;
        SS = 0;
        date(ii) = datenum(year(ii),month(ii),day(ii),hr(ii),MM,SS);
        twl_output(ii).tide = twl_output(ii).twl-twl_output(ii).runup;
    end

