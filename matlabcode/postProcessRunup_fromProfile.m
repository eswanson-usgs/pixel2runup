%Code to convert horizontal runup to elevation using topo profile

clear,close all
%--------------------------------------------------------------------------

%set paths (M.Palmsten 6/5/2021)
addpath \\gs\stpetersburgfl-g\NACCH\Code\matlab\m-cmg_backup\jbrown\TS_ANALYSIS
% addpath C:\Users\mpalmsten\Documents\MATLAB\trunk\stpete\noaa
% addpath C:\Users\mpalmsten\Documents\MATLAB\jbrown\RUNUP
%addpath('C:\Users\jbirchler\OneDrive - DOI\Documents\stpete_branches\jbrown\RUNUP')
addpath('\\gs\stpetersburgfl-g\NACCH\Code\matlab\m-cmg_backup\jbrown\TS_ANALYSIS')
addpath('\\gs\stpetersburgfl-g\NACCH\Code\matlab\stats')
addpath('\\gs\stpetersburgfl-g\NACCH\Code\matlab\lidar')
addpath('\\gs\stpetersburgfl-g\NACCH\Data\Molokai')
addpath('\\gs\stpetersburgfl-g\NACCH\Data\Molokai\code')
addpath('\\gs\stpetersburgfl-g\NACCH\Data\Molokai\code\ArgusTools\CILMatlab\nrlArgus\')
addpath('\\gs\stpetersburgfl-g\NACCH\Data\Molokai\code\ArgusTools\CILMatlab\pixel')
addpath('\\gs\stpetersburgfl-g\NACCH\Data\Molokai\code\ArgusTools\CILMatlab\argusDB')
addpath('\\gs\stpetersburgfl-g\NACCH\Data\Molokai\code\CIRN\CIL\CILTools\')
 %addpath('C:\Users\emilyjohnson\OneDrive - DOI\Documents\MATLAB\trunk')
 %addpath('C:\Users\emilyjohnson\OneDrive - DOI\Documents\MATLAB\trunk\stpete\noaa')
%addpath \\gs\stpetersburgfl-g\NACCH\Code\matlab\noaa

%***Inputs:
%---runup data:
runupFile.dates = [2018 08 05;
                   2018 08 16];
%---topo data:all
%-dates: 20180625;                         
topoFile.date = [2018 06 25]; 

topoFile.lineno = [03]; %this could also have an option to be an average, if 2 numbers are input it would use the mean profile
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
dates(:,1) = datenum(runupFile.dates(1,:)):datenum([0,0,1,0,0,0]):datenum(runupFile.dates(2,:));

% load database structure - in the future edit this to a call to the DB
load('\\gs\stpetersburgfl-g\NACCH\Data\Molokai\geometry\molokai.coastCam.mat');
cameraInfo = molokai;
clear molokai

for ii = 1:size(dates,1) %loop thru multiple days at once
    
runupFile.date = datevec(dates(ii));
runupFile.path = ['\\gs\stpetersburgfl-g\NACCH\Data\Molokai\runup\' num2str(runupFile.date(1))]; % location of digtized runup 
% runupFile.path = ['\\gs\stpetersburgfl-g\NACCH\Data\Molokai\runup\2018_gpsAlongshoreUniform\']; % location of digtized runup 
runupFile.stackPath = ['\\gs\stpetersburgfl-g\NACCH\Data\Molokai\argusArchive\', num2str(runupFile.date(1)), '\cx\' ]; % location of time stack

% topoFile.path = ['\\gs\stpetersburgfl-g\NACCH\Data\Molokai\geometry\' num2str(topoFile.date(1)) '' num2str(topoFile.date(2),'%2.2d') '' num2str(topoFile.date(3),'%2.2d')];
%topoFile.path = ['\\gs\stpetersburgfl-g\NACCH\Data\Molokai\geometry\'];
%[topoFile] = topoProfile(topoFile);

% read in interpolated survey structure surveyInterp created by
% importSurveyData.m
load(['\\gs\stpetersburgfl-g\NACCH\Data\Molokai\geometry\','topoSurvey.mat' ]);
surveyInterp.path = ['\\gs\stpetersburgfl-g\NACCH\Data\Molokai\geometry\','topoSurvey.mat' ];
%computeRunupElevation_fromProfile_fromTimestacks(runupFile,topoFile, surveyInterp);
computeRunupElevation_fromProfile_fromTimestacks(runupFile, cameraInfo, surveyInterp);

disp('All Done')
end
