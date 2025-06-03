%% load all runup data (timestacks)
% load 'runup' files, put them in structure
% now easily useable to compare to model TWL forecasts
clear;close all;clc

% load all_runup.mat file, load file directories
load('\\gs\StPetersburgFL-G\NACCH\Imagery\madbeach\runup\all_runupNew.mat')

count=length(all_runup);

year = 2021;
cd(['\\gs\StPetersburgFL-G\NACCH\Imagery\madbeach\runup\',num2str(year),'\'])
files = dir();
files = files(3:end);
files = files([files.isdir]);

% find last entry in all_runup, exclude runup already in file
tmp=(all_runup(end).name(1:8));
yrday=num2str(date_to_yearday(str2num(tmp(1:4)),str2num(tmp(5:6)),str2num(tmp(7:8))),'%03d');
monthstring=datestr([str2num(tmp(1:4)),str2num(tmp(5:6)),str2num(tmp(7:8)),0,0,0],'mmm');
day=(all_runup(end).name(7:8));
for ii=1:length(files)
    if strcmp(files(ii).name,[yrday,'_',monthstring,'.',day])
        [num2str(year),' - ',yrday,'_',monthstring,'.',day]
        files(1:ii)=[];
        break
    end
end

%%
for jj = 1:length(files)
    cd(['\\gs\StPetersburgFL-G\NACCH\Imagery\madbeach\runup\',num2str(year),'\',files(jj).name])
    yearday=str2num(files(jj).name(1:3));
    datev = datevec(datenum(year-1,12,31)+yearday);
    year = datev(1);
    month = datev(2);
    day = datev(3);
    
    runupfiles = dir([num2str(year),num2str(month,'%02g'),num2str(day,'%02g'),'*fromProfile.mat']);
    
    if ~isempty(runupfiles)
        for ii=1:length(runupfiles)
            runupfiles(ii).name
            load(runupfiles(ii).name)
            count=count+1;
            all_runup(count).name  = runupfiles(ii).name;
            all_runup(count).t     = Runup.t;
            all_runup(count).TWL   = Runup.TWL;
            all_runup(count).R2    = Runup.R2;
            all_runup(count).R5    = Runup.R5;
            all_runup(count).R10   = Runup.R10;
            all_runup(count).R13   = Runup.R13;
            all_runup(count).setup = Runup.setup;
            all_runup(count).S     = Runup.S;
            all_runup(count).Sinc  = Runup.Sinc;
            all_runup(count).SIG   = Runup.SIG;
            all_runup(count).slope = Runup.slope;
            all_runup(count).tide  = Runup.tide;
            
        end
    end
end

save('\\gs\StPetersburgFL-G\NACCH\Imagery\madbeach\runup\all_runupNew.mat','all_runup')


%%
if(0) % remove runup files based off OLD elevation profile
    for ii=1:length(all_runup)
        if all_runup(ii).t==all_runup(ii+1).t
            all_runup(ii)=[];
        end
    end
end

figure
plot(vertcat(all_runup.t),vertcat(all_runup.TWL),'x')
hold on
plot(vertcat(all_runup.t),vertcat(all_runup.R2),'o')
plot(vertcat(all_runup.t),vertcat(all_runup.tide),'s')
datetick('x','yyyy-mm-dd')
legend('TWL','R2','tide')
xlabel('date(yyyy-mm-dd)');ylabel('elevation')
title('Measured Water Levels from MadBeach Camera')


