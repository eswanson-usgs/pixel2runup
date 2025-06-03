%Convert horizontal runup to runup elevation
%   and computes R2% and setup
%   for runup timeseries from Argus timestacks
%
%Input:
%   runupName = runup timeseries filename
%
%   profile    = bathy/topo profile data (structure format)
%                   >> profile =
%                          x: [449881x1 double]
%                          y: [449881x1 double]
%                          z: [449881x1 double]
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
%   [Runup]=convertRunup_ArgusTimestack_madbeach(runupName,profile)
%--------------------------------------------------------------------------
function [Runup]=convertRunup_fromProfile_ArgusTimestack_madbeach(runupName,profile)

%--- madbeach
% inputFile = '\\gs\StPetersburgFL-G\NACCH\Imagery\madbeach\process\inputFile_20170213.m';
% eval(['run ',inputFile,';']); % this spits out the whole m-file...annoying
% geom = load(inputs.calib.extrinsic.geomFile);
% camera.x = geom.betas(1);
% camera.y = geom.betas(2);
% camera.z = geom.betas(3);


runupFile=parseFilename(runupName,'noLocal');

yloc    = str2double(runupFile.type(2:end)); %alongshore location of pixel transect
if isstring(runupFile.time) || ischar(runupFile.time)
    [year,mnth,day] = datevec(epoch2Matlab(str2double(runupFile.time)));
else
    [year,mnth,day] = datevec(epoch2Matlab(double(runupFile.time)));
end

monthName = datestr([year,mnth,day,0,0,0],'mmm');                              %3-letter abreviation of the month
dayStr    = num2str(day,'%2.2d');
yearday   = num2str(datenum(year,mnth,day)-datenum(year-1,12,31),'%3.3d');     %compute yearday (day of year, 0-365)


%---load runup timeseries (from Argus timestacks, created by extractRunupTimeseries_fromArgusTimestacks.m)
load (['\\gs\StPetersburgFL-G\NACCH\Imagery\madbeach\runup\' num2str(year) '\' yearday '_' monthName '.' dayStr '\' runupName])        %runup = {Ri; params; epoch; xyz; UV}

tRunup = runup.epoch-runup.epoch(1);    %time series (seconds)
etimeRunup = runup.epoch;               %epoch time
mtimeRunup = epoch2datenum(etimeRunup);  %matlab time (datenum)

iRunup = runup.Ri;                      %pixel indices of runup edge
xRunup = runup.xyz(:,1);                %x-coordinates of timestack pixels
yRunup = runup.xyz(:,2);                %y-coordinates of timestack pixels (constant)
zRunup = runup.xyz(:,3);


%---convert horizontal runup (x-dimension) to vertical runup (z-dimension) using measured topo data
%[xiRunup,yiRunup,ziRunup]=computeRunupMapped2Topo(xRunup,yRunup,zRunup,survey.x,survey.y,survey.z,camera.x,camera.y,camera.z);
ziRunup = interp1(profile.x,profile.z,xRunup);


%---Output:
Runup.yloc=yloc;
Runup.t=mean(mtimeRunup);

%***timeseries
Runup.ts=mtimeRunup;
Runup.x=xRunup;
Runup.y=yRunup;
Runup.xi=xRunup;
Runup.yi=yRunup;
Runup.zi=ziRunup;
