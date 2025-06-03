%% load all times camera was operating data
% load 'cx.runup90' files, put them in structure
% now easily useable to compare to model TWL forecasts
clear; close all; clc

% load all_camop.mat file, load file directories
load('\\gs\StPetersburgFL-G\NACCH\Imagery\madbeach\runup\all_camop.mat')
count=length(all_camop);
year = 2022;
cd(['\\gs\stpetersburgfl-g\Coastal_Change_Hazards\Archive\Data\2016\2016-363-DD_20161028\madbeach\',num2str(year),'\cx\'])
files = dir();
files = files(3:end);
files = files([files.isdir]);

for ii=1:length(files)
    k=strcmp(files(ii).name,all_camop(end).name);
    if k==1
        kk=ii
        return
    else
        kk=1;
    end
end

for ii=kk:length(files)
    jj=ii-kk+count;
all_camop(jj).name = files(ii).name;
all_camop(jj).folder = files(ii).folder;
day = str2num(all_camop(jj).name(9:10));
m = month(datetime([(all_camop(jj).name(9:10)),'-',all_camop(jj).name(5:7),'-',num2str(year)]));
all_camop(jj).date = datenum([year,m,day]);
end

% might need to remove some rows in case it captures some dates twice

save('\\gs\StPetersburgFL-G\NACCH\Imagery\madbeach\runup\all_camop.mat','all_camop')

% 2017 -> 343 items -> 343/345 days -> missing 2  (cumsum 2)
% 2018 -> 353 items -> 353/365 days -> missing 12 (cumsum 14)
% 2019 -> 356 items -> 356/365 days -> missing 9  (cumsum 23)
% 2020 -> 358 items -> 358/366 days -> missing 8  (cumsum 31)
% 2021 -> 364 items -> 364/365 days -> missing 1  (cumsum 32)
% 2022 -> items -> /365 days -> missing 
days_tot = round(datenum(clock()))-datenum(2017,1,21);
days_meas = length(all_camop);
days_miss = days_tot - days_meas;

figure(1);clf
vline(datenum(2017,1,21):round(datenum(clock())),'r')
hold on; box on
vline(vertcat(all_camop(:).date),'g')
title('Madeira Beach Camera Data Availability')
xlabel('Date')
datetick('x','yyyy-mm-dd')
set(gca,'ytick',[])
axis([datenum(2017,1,16) datenum(clock())+5 -0.01 1.01])
% set(gca,'xtick',[datenum(2017,1,1) datenum(2017,6,1) datenum(2018,1,1) datenum(2018,6,1) datenum(2019,1,1) datenum(2019,6,1) datenum(2020,1,1) datenum(2020,6,1) datenum(2021,1,1)])
text(datenum(clock())-210,.80,['Data start = ',datestr(all_camop(1).date)])
text(datenum(clock())-210,.76,['Data end  = ',datestr(all_camop(end).date)])
text(datenum(clock())-210,.725,['Total Days           = ',num2str(days_tot)])
text(datenum(clock())-210,.69, ['Days with data     = ',num2str(days_meas)])
text(datenum(clock())-210,.65, ['Days with no data = ',num2str(days_miss)])


%% find last entry in all_camop, exclude runup already in file
tmp=datevec((all_camop(end).date));
% yrday=num2str(date_to_yearday(str2num(tmp(1:4)),str2num(tmp(5:6)),str2num(tmp(7:8))),'%03d');
% monthstring=datestr([str2num(tmp(1:4)),str2num(tmp(5:6)),str2num(tmp(7:8)),0,0,0],'mmm');
yrday=num2str(date_to_yearday((tmp(1)),(tmp(2)),(tmp(3))));
monthstring=datestr([(tmp(1)),(tmp(2)),(tmp(3)),0,0,0],'mmm');
day=(all_camop(end).name(9:10));
for ii=1:length(files)
    if strcmp(files(ii).name,[yrday,'_',monthstring,'.',day])
        [yrday,'_',monthstring,'.',day]
        files(1:ii)=[];
        return
    end
end

%% 
for jj = 1:length(files)
    count=count+1;
    all_camop(count).name  = files(jj).name;
    all_camop(count).folder  = files(jj).folder;
    all_camop(count).date  = datenum([(files(jj).folder(end-6:end-3)),...
        (files(jj).name(end-5:end-3)),(files(jj).name(end-1:end))],'yyyymmmdd');
end

% for jj = 1:length(files)
%     cd(['\\gs\StPetersburgFL-G\NACCH\Imagery\madbeach\runup\',num2str(year),'\',files(jj).name])
%     yearday=str2num(files(jj).name(1:3));
%     datev = datevec(datenum(year-1,12,31)+yearday);
%     year = datev(1);
%     month = datev(2);
%     day = datev(3);
%     
%     runupfiles = dir([num2str(year),num2str(month,'%02g'),num2str(day,'%02g'),'*fromProfile.mat'])
%     
%     if ~isempty(runupfiles)
%         for ii=1:length(runupfiles)
%             runupfiles(ii).name
%             load(runupfiles(ii).name)
%             count=count+1;
%             all_camop(count).name  = runupfiles(ii).name;
%             all_camop(count).t     = Runup.t;
%             all_camop(count).TWL   = Runup.TWL;
%             all_camop(count).R2    = Runup.R2;
%             all_camop(count).R5    = Runup.R5;
%             all_camop(count).R10   = Runup.R10;
%             all_camop(count).R13   = Runup.R13;
%             all_camop(count).setup = Runup.setup;
%             all_camop(count).S     = Runup.S;
%             all_camop(count).Sinc  = Runup.Sinc;
%             all_camop(count).SIG   = Runup.SIG;
%             all_camop(count).slope = Runup.slope;
%             all_camop(count).tide  = Runup.tide;
%             
%         end
%     end
% end

save('\\gs\StPetersburgFL-G\NACCH\Imagery\madbeach\runup\all_camop.mat','all_camop')


