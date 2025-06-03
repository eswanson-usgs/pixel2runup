%% 
%Extract runup from madbeach timestacks using 'runupTool_madbeach.m'
%
%Input:
%   dataDIR = full directory path where timestacks are located
%   saveDIR = full directory path where timeseries are saved
%
%   - date info
%   - alongshore locations of runup timestacks
%
%Output:
%   - saves runup timeseries data to .mat file (same name as timestack with .runup)
%
%--------------------------------------------------------------------------
clear,close all
clc

%***Input:
year=2022;
month= 09;
day  = 09;
yrunup=[90]; %[90, 150, 250] = alongshore locations of runup pixel lines/timestacks

%---timestack files located at:
% dataDIR = ['\\gs\stpetersburgfl-g\Coastal_Change_Hazards\Archive\Data\2022\2022-305-FA\madbeach\',num2str(year),'\cx'];
dataDIR = ['\\gs\stpetersburgfl-g\Coastal_Change_Hazards\Archive\Data\2022\2022-305-FA\madbeach\',num2str(year),'\c1'];
%---save runup timeseries files to:
saveDIR = ['\\gs\StPetersburgFL-G\NACCH\Imagery\madbeach\runup\',num2str(year)];%,'\'];

%--------------------------------------------------------------------------

yearday=num2str(datenum(year,month,day)-datenum(year-1,12,31),'%3.3d');
monthName=datestr([year,month,day,0,0,0],'mmm');
dayStr=num2str(day,'%2.2d');
dirName = [yearday '_' monthName '.' dayStr];

dataPath = [dataDIR '\' dirName ];
savePath = [saveDIR '\' dirName];
if ~exist(savePath)
    mkdir(savePath)
end

%---for each alongshore runup location:
for yy=1:length(yrunup)
    
    S = dironly([dataPath '\*runup' num2str(yrunup(yy),'%03d') '.mat']);
    
    for ss=1:length(S)
        
        figure('Name',S(ss).name);
        S(ss).name
        runupTool_madbeach([dataPath '\' S(ss).name],savePath)
        
        if ss==length(S)
            keyboard
        end
    end
end

%% REDO RUNUP EXTRACTION TO CORRECT FOR ERRORS
 if (0)
      
     %      set(0, 'DefaultFigurePosition', [50 50 1000 600]);
%      set(0, 'DefaultFigurePosition', [1375 -350 1250 1000]); %jjb telework default (2020,laptop/big screen)
     set(0,'DefaultFigurePosition',[1925 -230 1250 1000]); %jjb new telework default (2021,laptop/big screen)
     
     clear,close all
     
%      cd \\gs\StPetersburgFL-G\NACCH\Imagery\madbeach\runup\2018\189_Jul.08\
     cd \\gs\StPetersburgFL-G\NACCH\Imagery\madbeach\runup\2017\280_Oct.07\
     
     
     files = dir('*madbeach.cx.*.runup.mat');
     savePath = pwd
     
     for  ii = 1:length(files)
         imageName = files(ii).name
         runupTool_madbeach(imageName,savePath)
         keyboard
     end
      
 end

%% REDO RUNUP - justin
if(0)

% %     set(0, 'DefaultFigurePosition', [1923 110 1200 1000]);
% %     set(0, 'DefaultFigurePosition', [4 44 1019 640]);
%     set(0,'DefaultFigurePosition',[1367 -56 1084 740]);
set(0,'DefaultFigurePosition',[1925 -295 1220 1075])
clear; close all

    cd \\gs\StPetersburgFL-G\NACCH\Imagery\madbeach\runup\2020
    folders = dir();
    for jj=192:349 % looking back on july2020

% cd \\gs\StPetersburgFL-G\NACCH\Imagery\madbeach\runup\2022
% folders = dir();
% for jj=23:-1:3 % 2022
    
    cd(folders(jj).name)
    %         files = dir('*madbeach.cx.runup90.runup.mat');
    files = dir('*madbeach.cx.r*90.runup.mat');
    savePath = pwd
    for  ii = 1:length(files)
        imageName = files(ii).name
        runupTool_madbeach(imageName,savePath)
        keyboard
    end
    cd ..
end
    
end
