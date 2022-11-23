%Code to convert horizontal runup to elevation using topo profile

clear,close all
%--------------------------------------------------------------------------
%***Inputs:

%---runup data:
% runupFile.dates = [2017 07 15;
%                    2017 07 15];
 runupFile.dates = [2020 09 16;
                    2020 09 16];


%---topo data:
%all dates: 20160909; 20161130; 20170217; 20170509; 20170914; 20171109; 20180124
%           20181015; 20190918; 20200610; 20200710; 20200908; 20200921; 20201106; 
%           20201116; 20201218; 20210115; 20210303; 20210421; 20210616;
%-dates DONT have line 60: 20160909; 20161130; 
%                          20170217; 20170914;
%-dates      HAVE line 60: 20170509; 20171109; 
%                          20180124; 20181015; 
%                          20190918; 
%                          20200610; 20200710; 20200908; 20200921; 20201106; 
%                          20201116; 20201218;
%                          20210115; 20210303; 20210421; 20210616;

%                                       
%topoFile.date = [2017 05 09]; 
topoFile.date = [2020 09 08]; 

% choose line 60 if available, or choose to average 48/49
topoFile.lineno = [60]; %this could also have an option to be an average, like if 2 numbers are input it would use the mean profile
% topoFile.lineno = [48 49]; %this takes average of lines 48/49, closest lines to 60
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
dates(:,1) = datenum(runupFile.dates(1,:)):datenum([0,0,1,0,0,0]):datenum(runupFile.dates(2,:));

for ii = 1:size(dates,1) %loop thru multiple days at once
    
runupFile.date = datevec(dates(ii));
runupFile.path = ['\\gs\StPetersburgFL-G\NACCH\Imagery\madbeach\runup\' num2str(runupFile.date(1))];

topoFile.path = ['\\gs\StPetersburgFL-G\NACCH\Imagery\madbeach\surveys\walking\' num2str(topoFile.date(1)) '' num2str(topoFile.date(2),'%2.2d') '' num2str(topoFile.date(3),'%2.2d')];

[topoFile] = topoProfile_madbeach_y90(topoFile);

computeRunupElevation_fromProfile_fromTimestacks_madbeach(runupFile,topoFile);

disp('All Done')
end