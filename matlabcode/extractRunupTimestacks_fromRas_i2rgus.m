%%
%Extract runup timestacks from madbeach i2rgus ras.tiff file
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
camera = 'c2'

%---timestack files located at:
dataDIR = ['\\gs\stpetersburgfl-g\Coastal_Change_Hazards\Archive\Data\2022\2022-305-FA\madbeach\',num2str(year),'\',camera];
%---save runup timeseries files to:
saveDIR = dataDIR;

yearday=num2str(datenum(year,month,day)-datenum(year-1,12,31),'%3.3d');
monthName=datestr([year,month,day,0,0,0],'mmm');
dayStr=num2str(day,'%2.2d');
dirName = [yearday '_' monthName '.' dayStr];

dataPath = [dataDIR '\' dirName ];
savePath = [saveDIR '\' dirName];


%% make sure you know how the ras.tiff file is built; i.e. which timestacks are part of file
% see 'madbeach\runup\pix' for pix files and README
clear y20 y90 y150 y250
if datenum(year,month,day)>datenum(2022,08,04) && datenum(year,month,day)<datenum(2022,09,09,13,30,0)
    if strcmp(camera,'c1')
        y20  =    1: 350;
        y90  =  351:1198;
        y150 = 1199:1988;
        y250 = 1988:2682;
    elseif strcmp(camera,'c2')
        y20  =    1: 622;
    end
load('\\gs\stpetersburgfl-g\NACCH\Imagery\madbeach\runup\madbeach-runup\old_UTM.XYZ\runup20_UTM.XYZ.mat')
load('\\gs\stpetersburgfl-g\NACCH\Imagery\madbeach\runup\madbeach-runup\old_UTM.XYZ\runup90_UTM.XYZ.mat')
load('\\gs\stpetersburgfl-g\NACCH\Imagery\madbeach\runup\madbeach-runup\old_UTM.XYZ\runup150_UTM.XYZ.mat')
load('\\gs\stpetersburgfl-g\NACCH\Imagery\madbeach\runup\madbeach-runup\old_UTM.XYZ\runup250_UTM.XYZ.mat')
elseif datenum(year,month,day)>datenum(2022,09,09,13,0,0)
    if strcmp(camera,'c1')
        y20  =    1: 699;
        y90  =  700:2394;
        y150 = 2395:3972;
        y250 = 3973:5358;
    elseif strcmp(camera,'c2')
        y20  =    1:1243;
    end
load('\\gs\stpetersburgfl-g\NACCH\Imagery\madbeach\runup\madbeach-runup\runup20_UTM.XYZ.mat')
load('\\gs\stpetersburgfl-g\NACCH\Imagery\madbeach\runup\madbeach-runup\runup90_UTM.XYZ.mat')
load('\\gs\stpetersburgfl-g\NACCH\Imagery\madbeach\runup\madbeach-runup\runup150_UTM.XYZ.mat')
load('\\gs\stpetersburgfl-g\NACCH\Imagery\madbeach\runup\madbeach-runup\runup250_UTM.XYZ.mat')
end

%% now load ras.tiff and split into individual files
cd(dataPath)
files_runup = dir(['*',camera,'.ras.tiff']);
for ii=1:length(files_runup)
    tiff_runup = imread(files_runup(ii).name);
    if exist('y20','var')
        RAW = tiff_runup(:,y20,:);
        T = str2num(files_runup(ii).name(1:10)):0.5:str2num(files_runup(ii).name(1:10))+((size(tiff_runup,1)-1)/2);
        XYZ = [runup20.XYZ(:,1),runup20.XYZ(:,2),zeros(length(runup20.XYZ),1)];
        save([files_runup(ii).name(1:end-8),'cx.runup020.mat'],'RAW','T','XYZ')
        clear RAW T XYZ
    end
    if exist('y90','var')
        RAW = tiff_runup(:,y90,:);
        T = str2num(files_runup(ii).name(1:10)):0.5:str2num(files_runup(ii).name(1:10))+((size(tiff_runup,1)-1)/2);
        XYZ = [runup90.XYZ(:,1),runup90.XYZ(:,2),zeros(length(runup90.XYZ),1)];
        save([files_runup(ii).name(1:end-8),'cx.runup090.mat'],'RAW','T','XYZ')
        clear RAW T XYZ
    end
    if exist('y150','var')
        RAW = tiff_runup(:,y150,:);
        T = str2num(files_runup(ii).name(1:10)):0.5:str2num(files_runup(ii).name(1:10))+((size(tiff_runup,1)-1)/2);
        XYZ = [runup150.XYZ(:,1),runup150.XYZ(:,2),zeros(length(runup150.XYZ),1)];
        save([files_runup(ii).name(1:end-8),'cx.runup150.mat'],'RAW','T','XYZ')
        clear RAW T XYZ
    end
    if exist('y250','var')
        RAW = tiff_runup(:,y250,:);
        T = str2num(files_runup(ii).name(1:10)):0.5:str2num(files_runup(ii).name(1:10))+((size(tiff_runup,1)-1)/2);
        XYZ = [runup250.XYZ(:,1),runup250.XYZ(:,2),zeros(length(runup250.XYZ),1)];
        save([files_runup(ii).name(1:end-8),'cx.runup250.mat'],'RAW','T','XYZ')
        clear RAW T XYZ
    end
end


