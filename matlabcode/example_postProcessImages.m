clear,close all
clc
%run C:\Users\jenniferbrown\Documents\MATLAB\SVN\stpete\nctoolbox\setup_nctoolbox.m
% output dir = \\gs\StPetersburgFL-G\NACCH\Imagery\madbeach\2017\data\
inputFile = '\\gs\StPetersburgFL-G\NACCH\Imagery\madbeach\process\inputFile_20170213.m';
ext='jpg';

%*** load inputFile
eval(['run ' inputFile ]);

%--------------------------------------------------------------------------
%--- 0. Date Info
% time = clock;
% year = time(1);
% month = time(2);
% day = time(3);
% yearday=datenum(year,month,day)-datenum(year-1,12,31);
% 
% %---get yesterday's date
% [yearday0,year0,month0,day0]=yesterday(year,yearday);
% yearday0 = num2str(yearday0,'%3.3d');
% 
% monthName=datestr([year0,month0,day0,0,0,0],'mmm');
% dayStr=num2str(day0,'%2.2d');
% dirName = [yearday0 '_' monthName '.' dayStr];

%***loop through multiple days
year = 2020;
month = 9;
days = [16:17];

for dd = 1:length(days)
    
yearday=datenum(year,month,days(dd))-datenum(year-1,12,31);
yearday=num2str(yearday,'%03d');
monthName=datestr([year,month,days(dd),0,0,0],'mmm');
dayStr=num2str(days(dd),'%2.2d');
dirName = [yearday '_' monthName '.' dayStr];
    
% %--------------------------------------------------------------------------
%--- 1. Image Products
dirIN = [inputs.path.c1 '\' dirName];
dirIN = ['\\gs\StPetersburgFL-G\NACCH\Archive\Data\2016\2016-363-DD_20161028\madbeach\' num2str(year) '\c1\' dirName];
dirOUT = [inputs.path.data '\' dirName];
dirOUT = ['\\gs\StPetersburgFL-G\NACCH\Imagery\madbeach\' num2str(year) '\data\' dirName];

postProcessImages_madbeach_20170213(dirIN,dirOUT,inputs,ext);

clear dirIN dirOUT
%--------------------------------------------------------------------------
%--- 2. Pixel Products

% dirIN = [inputs.path.cx '\' dirName];
% dirOUT = [inputs.path.data '\' dirName];
% 
% postProcessPixelProducts_madbeach_20170213(dirIN,dirOUT,inputs,ext);
% 
% clear dirIN dirOUT
%--------------------------------------------------------------------------
clear yearday monthName dayStr dirName
end
