%Convert horizontal runup to runup elevation
%   and computes R2% and setup
%   for runup timeseries from Argus timestacks
%
%Input:
%   runupName = runup timeseries filename
%
%   survey    = bathy/topo profile data (structure format)
%                   >> profile =
%                          x: [449881x1 double]
%                          y: [449881x1 double]
%                          z: [449881x1 double]
%   cameraInfo = structure of camera information from coastCam data base
%
%Output:
%   Runup = runup timeseries data (structure format)
%               >> Runup = 
%                      yloc: 1100
%                         t: 7.3604e+05
%                        ts: [2048x1 double]
%                         x: [2048x1 double]
%                         y: [2048x1 double]
%                        xi: [2048x1 double]
%                        yi: [2048x1 double]
%                        zi: [2048x1 double]
%
%Usage: 
%   [Runup]=convertRunup_ArgusTimestack(runupName,profile)
%--------------------------------------------------------------------------
function [Runup]=convertRunup_fromProfile_ArgusTimestack(runupName,survey,cameraInfo)

%--- madbeach
% inputFile = 'C:\Imagery\process\inputFile_20170213.m';
% eval(['run ',inputFile,';']); % this spits out the whole m-file...annoying
% geom = load(inputs.calib.extrinsic.geomFile);

runupFile=parseFilename(runupName.name,'noLocal');

yloc    = str2double(runupFile.type(end-1:end)); %alongshore location of pixel transect
if isstring(runupFile.time) || ischar(runupFile.time)
    [year,mnth,day] = datevec(epoch2Matlab(str2double(runupFile.time)));
else
    [year,mnth,day] = datevec(epoch2Matlab(double(runupFile.time)));
end

monthName = datestr([year,mnth,day,0,0,0],'mmm');                              %3-letter abreviation of the month
dayStr    = num2str(day,'%2.2d');
yearday   = num2str(datenum(year,mnth,day)-datenum(year-1,12,31),'%3.3d');     %compute yearday (day of year, 0-365)

%---load runup timeseries (from Argus timestacks, created by extractRunupTimeseries_fromArgusTimestacks.m)
load ([runupName.path '\' yearday '_' monthName '.' dayStr '\' runupName.name])        %runup = {Ri; params; epoch; xyz; UV}

% load runup stack because it contains (x, y, 0) values of stack
runupName.stackName = [runupName.time, '.', runupName.when, '.',runupName.station,'.c',runupName.camera, '.', runupName.type, '.mat'];
load ([runupName.stackPath  '\' yearday '_' monthName '.' dayStr '\' runupName.stackName])

tRunup = runup.epoch-runup.epoch(1);    %time series (seconds)
etimeRunup = runup.epoch;               %epoch time
mtimeRunup = epoch2datenum(etimeRunup);  %matlab time (datenum)

% (x,y, 0) location of digtized runup
iRunup = runup.Ri;                      %pixel indices of runup edge
xRunup = XYZ(runup.Ri,1);
yRunup = XYZ(runup.Ri,2);
zRunup = XYZ(runup.Ri,3);
%
%---convert horizontal runup (x-dimension) to vertical runup (z-dimension) using measured topo data
[xiRunup,yiRunup,ziRunup]=computeRunupMapped2Topo(xRunup,yRunup,zRunup,survey.x,survey.y,survey.z,cameraInfo.camera.x,cameraInfo.camera.y,cameraInfo.camera.z);

%---Output:
Runup.yloc=yloc;
Runup.t=mean(mtimeRunup);

%***timeseries
Runup.ts=mtimeRunup;
Runup.x=xRunup;
Runup.y=yRunup;
Runup.xi=xiRunup;
Runup.yi=yiRunup;
Runup.zi=ziRunup;
