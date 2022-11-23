%% load past TWL forecasts
% 
clear; close all; clc
MySQLDBConnect%('IGSAFPESGSZ03','3306','twlViewer','nacch','twl763')
%tampa bay region
region = 'TBW';
% position of camera
% latitude  =  27.796216949206798; 
% longitude = -82.796102542950635;
latitude    =  27.7962; % MySQL DB having trouble with too many decimals
longitude   = -82.7961;
%date
dateTime='2017-01-21 00:00:00';
time1=datenum(dateTime);

% load in previous data, set up to save new data to file
load('\\gs\StPetersburgFL-G\NACCH\Imagery\madbeach\runup\all_TWL_forecast.mat')
% count=length(all_twl);
year = 2017;
cd(['\\gs\StPetersburgFL-G\NACCH\Imagery\madbeach\runup\'])

timeAll = datenum(vertcat(all_twl(:).forecastTime));
time1=datenum(all_twl(end).forecastTime);
dateTime = datestr(time1,'yyyy-mm-dd HH:MM:SS')
ii=0

%%
for ii=1:1000

    time1 = time1+(1/24);
    dateTime = datestr(time1,'yyyy-mm-dd HH:MM:SS');
    
    if ceil(time1)-time1<(1/240)
        ii
        dateTime
        time1=round(time1);
    end
    
    if time1 > datenum(clock)-1
        display('dont load any more data: last 24hrs have been reached')
        display(['saving'])
        save('\\gs\StPetersburgFL-G\NACCH\Imagery\madbeach\runup\all_TWL_forecast.mat','all_twl')
        display(' ')
        display(['done-final dateTime = ',all_twl(end).forecastTime])
        return
    end
    
    % this will load all forecasts within the ~104 hr window before/after dateTime
    % [twl_output] = MySQLDBGetTWLTimeseries(longitude, latitude, dateTime, region);
    
    % this will load all forecasts for dateTime ONLY (12 timepoints?)
    [twl_output] = MySQLDBGetTWLForecast(longitude, latitude, dateTime, region);
    if ~isfield(twl_output,'slope')
        twl_output(end).slope=NaN;
    end
    if ~isfield(twl_output,'hs')
        twl_output(end).hs=NaN;
    end
    if ~isfield(twl_output,'tp')
        % period as of ~5/2020 is labeled 'pp' not 'tp'
        [twl_output.tp] = twl_output.pp; 
        twl_output = orderfields(twl_output,[1:21,28,22:27]); 
        twl_output = rmfield(twl_output,'pp');
%         twl_output(end).tp=NaN;
    end
%     if ~isfield(twl_output,'tide')
%         twl_output(end).tide=NaN;
%     end
%     if ~isfield(twl_output,'surge')
%         twl_output(end).surge=NaN;
%     end
    all_twl = [all_twl;twl_output(end)];
    % keyboard
    if ii==300 || ii==600 %|| ii==600 || ii==700 || ii==900
        display(['saving: thru ',num2str(ii),' timesteps'])
        save('\\gs\StPetersburgFL-G\NACCH\Imagery\madbeach\runup\all_TWL_forecast.mat','all_twl')
    end
end
% if ii/500==round(ii/500)
display(' ')
display(['done-saving'])
save('\\gs\StPetersburgFL-G\NACCH\Imagery\madbeach\runup\all_TWL_forecast.mat','all_twl')
display(' ')
display(['done-final dateTime = ',all_twl(end).forecastTime])
% end
% for ii = 1:length(twl_output)
% year(ii) = str2num(twl_output(ii).forecastTime(1:4));
% month(ii) = str2num(twl_output(ii).forecastTime(6:7));
% day(ii) = str2num(twl_output(ii).forecastTime(9:10));
% hr(ii) = str2num(twl_output(ii).forecastTime(12:13));
% MM = 0;
% SS = 0;
% date(ii) = datenum(year(ii),month(ii),day(ii),hr(ii),MM,SS);
% twl_output(ii).tide = twl_output(ii).twl-twl_output(ii).runup;
% end
% % only use this if wanting to compare "forecastBegins" times:
% for ii = 1:length(twl_output)
% year(ii) = str2num(twl_output(ii).forecastBegins(1:4));
% month(ii) = str2num(twl_output(ii).forecastBegins(6:7));
% day(ii) = str2num(twl_output(ii).forecastBegins(9:10));
% hr(ii) = str2num(twl_output(ii).forecastBegins(12:13));
% MM = 0;
% SS = 0;
% date(ii) = datenum(year(ii),month(ii),day(ii),hr(ii),MM,SS);
% twl_output(ii).tide = twl_output(ii).twl-twl_output(ii).runup;
% end

