%% Plot TWL estimates on raw camera image
% also plots measured runup values on camera image
% plots cross-shore profile with wl values along it
% cd \\gs\StPetersburgFL-G\NACCH\Imagery\madbeach\
% clear,close all
% addpath(genpath('C:\Users\jwlong\Documents\MATLAB\m-cmg\branches\jbrown\CIRN\CIL\UAV-Processing-Toolbox'))

%add paths
addpath \\gs\stpetersburgfl-g\NACCH\Data\TresPalmasPR\code
addpath \\gs\stpetersburgfl-g\NACCH\Imagery\madbeach\matlab
addpath \\gs\stpetersburgfl-g\NACCH\Imagery\madbeach\process
addpath \\gs\stpetersburgfl-g\NACCH\Imagery\madbeach\process\runupConvertToVertical

if strcmp('jbirchler',getenv('username'))
addpath(genpath('C:\Users\jbirchler\OneDrive - DOI\Documents\stpete_branches\jbrown\CIRN\CIL\UAV-Processing-Toolbox'))
end


year  = 2022;
month = 09;  monthNum=num2str(month,'%02d');
day   = 15;  dayNum=num2str(day,'%02d');  
hr1   = 13;
% hr2  = 18;
yearday=datenum(year,month,day)-datenum(year-1,12,31);yearday=num2str(yearday,'%03d');
monthName=datestr([year,month,day,0,0,0],'mmm');
myepoch = num2str(datenum2epoch(datenum(year,month,day,hr1,0,0)));
[~,dayName] = weekday(datenum(year,month,day));


%*** 1. image data
%load(['\\gs\StPetersburgFL-G\Coastal_Change_Hazards\Archive\Data\2016\2016-363-DD_20161028\madbeach\',num2str(year),'\c1\',yearday,'_',monthName,'.',dayNum,'\',myepoch,'.',dayName,'.',monthName,'.',dayNum,'_',num2str(hr1),'_00_00.GMT.',num2str(year),'.madbeach.c1.snap.mat']);
%load(['\\gs\StPetersburgFL-G\Coastal_Change_Hazards\Archive\Data\2016\2016-363-DD_20161028\madbeach\',num2str(year),'\c1\',yearday,'_',monthName,'.',dayNum,'\',myepoch,'.',dayName,'.',monthName,'.',dayNum,'_',num2str(hr1),'_00_00.GMT.',num2str(year),'.madbeach.c1.timex.mat']);
geom = load('geomFile_c1.mat');

figure(1);
image1 = imread('1663246800.Thu.Sep.15_13_00_00.GMT.2022.madbeach.c1.snap.jpg');
imagesc(image1)
%imshow('1663246800.Thu.Sep.15_13_00_00.GMT.2022.madbeach.c1.snap.jpg');
%print('-dpng',['C:\Users\mpalmsten\OneDrive - DOI\projects\TWL\MadeiraBeach\ForecastOnImage\Snap',num2str(year), num2str(month), num2str(day), num2str(hr1), '.png'])

figure(2);
image2 = imread('1663246800.Thu.Sep.15_13_00_00.GMT.2022.madbeach.c1.timex.jpg');
imagesc(image2)
%imshow('1663246800.Thu.Sep.15_13_00_00.GMT.2022.madbeach.c1.timex.jpg')
% print('-dpng', ['C:\Users\mpalmsten\OneDrive - DOI\projects\TWL\MadeiraBeach\ForecastOnImage\Timex',num2str(year), num2str(month), num2str(day), num2str(hr1), '.png'])
% 
% *** 2. topo data
% load('\\gs\StPetersburgFL-G\NACCH\Imagery\madbeach\surveys\walking\20170217\20170217Local.mat','topo');
% load('\\gs\StPetersburgFL-G\NACCH\Imagery\madbeach\surveys\walking\20170509\20170509_topoProfile_y90.mat');
% surveys = ['20160909';'20161130';'20170217';'20170509';'20170914';'20171109';'20180124'];
surveys_all = dir('\\gs\stpetersburgfl-g\NACCH\Imagery\madbeach\surveys\walking\2*');
surveys = [];
for ii=1:length(surveys_all); surveys = [surveys; surveys_all(ii).name];end %specify
tmp=((datenum(year,month,day)-datenum(surveys,'yyyymmdd')));
ii=find(tmp>0,1,'last'); %find most recent topo
if ii==4 || ii>=6
ans = dlmread(['\\gs\StPetersburgFL-G\NACCH\Imagery\madbeach\surveys\walking\20221005\line60.xyz']); %%%%%%%changing to 9/15/22 for testing
profile.E=ans(:,1);profile.N=ans(:,2);
[profile.x,profile.y]=coordSys_madbeach(profile.E,profile.N);
profile.z=ans(:,3);
elseif ii==1 || ii==2 || ii==3 || ii==5
disp('hi')
dlmread(['\\gs\StPetersburgFL-G\NACCH\Imagery\madbeach\surveys\walking\20221005\line60.xyz']);
p1.E=ans(:,1);p1.N=ans(:,2);
[p1.x,p1.y]=coordSys_madbeach(p1.E,p1.N);
p1.z=ans(:,3);
dlmread(['\\gs\StPetersburgFL-G\NACCH\Imagery\madbeach\surveys\walking\20221005\line60.xyz']);
p2.E=ans(:,1);p2.N=ans(:,2);
[p2.x,p2.y]=coordSys_madbeach(p2.E,p2.N);
p2.z=ans(:,3);
p2.x_new=p1.x;
p2.z_new= interp1(p2.x,p2.z,p2.x_new);
profile.x=p1.x;  profile.z=mean([p1.z, p2.z_new],2);
profile.y=ones(size(p1.y))*-87.7; % estimating line 60 position to be -87.7 based on surveys that measured line 60
end

profile.y=ones(size(profile.y))*-87.7; % estimating line 60 position to be -87.7 based on surveys that measured line 60

profile.z(isnan(profile.z))=0; %change NaNs to 0s
[a,b]=histcounts(profile.z,unique(profile.z));
mult=find(a>1);
notUnique = (size(profile.z) ~= size(unique(profile.z)));
while notUnique(1) == 1 % wasnt catching >2 identical values
    for i=1:length(mult)
        %b is ascending array of unique values. mult is array of bins in a with dupllicate values. indices in a match indices in b. So, index mult(i) will give index of a value in b that is not unique
        dupValue = b(mult(i));
        dupIndices = find(profile.z == dupValue);

        for ii=1:length(dupIndices)
            %dupIndices(ii) is index in profile.z that contains duplicate values
            index = dupIndices(ii);
            %slightly alter value to prevent duplicates
            offset = rand * 0.0000000001;
            profile.z(index) = profile.z(index) + offset;
        end
    end
    [a,b]=histcounts(profile.z,unique(profile.z));
    mult=find(a>1);
    notUnique = (size(profile.z) ~= size(unique(profile.z)));
end

% % yindex = find(min(abs(topo.y(:,1)-(-90)))==abs(topo.y(:,1)-(-90)));
% % profile.x = topo.x(yindex,:);
% % profile.z = topo.z(yindex,:);
% % good = find(isnan(profile.z)==0);
% % profile.x = profile.x(good);
% % profile.z = profile.z(good);
% % profile.y = -90.*ones(size(profile.x));
% 
% % figure(10);
% % plot(profile.x,profile.z,'k');
% 
% 
% %*** 3. model data
cam.lat = 27.796216949206798;
cam.lon = -82.796102542950635;
% % load(['\\gs\StPetersburgFL-G\NACCH\Projects\NWS\OutputArchive\WCOSS_test\monthFiles\',num2str(month),'_2017runup.mat']);
% % pindex = find(min(abs(R.latitude-cam.lat))==abs(R.latitude-cam.lat));
% % % pindex = find(min(abs(R.Lat_Update-cam.lat))==abs(R.Lat_Update-cam.lat));
% % tindex = find(R.time>datenum(year,month,day,hr1,0,0) & R.time<datenum(year,month,day,hr2,0,0));
% 
load('\\gs\StPetersburgFL-G\NACCH\Imagery\madbeach\runup\all_TWL_forecast.mat')
R.time=datenum(cat(1,all_twl.forecastTime),'yyyy-mm-dd HH:MM:SS.FFF')';
R.runup05=cat(1,all_twl.runup05)';
R.runup=cat(1,all_twl.runup)';
R.runup95=cat(1,all_twl.runup95)';
R.twl05=cat(1,all_twl.twl05)';    
R.twl=cat(1,all_twl.twl)';    
R.twl95=cat(1,all_twl.twl95)';
R.slope=cat(1,all_twl.slope)';

pindex = 1; % this file only saves the MadBeach TWL data
tindex = find(R.time>datenum(year,month,day,hr1-1,0,0) & R.time<datenum(year,month,day,hr1+1,0,0));

%Aadded by Eric - 11/23. Manually set instrincs array from meta.globals.lcp
%struct as input to xyzDistUV()
intrinsics = zeros(1,11);
intrinsics(1) = meta.globals.lcp.NU;
intrinsics(2) = meta.globals.lcp.NV;
intrinsics(3) = meta.globals.lcp.c0U;
intrinsics(4) = meta.globals.lcp.c0V;
intrinsics(5) = meta.globals.lcp.fx;
intrinsics(6) = meta.globals.lcp.fy;
intrinsics(7) = meta.globals.lcp.d1;
intrinsics(8) = meta.globals.lcp.d2;
intrinsics(9) = meta.globals.lcp.d3;
intrinsics(10) = meta.globals.lcp.t1;
intrinsics(11) = meta.globals.lcp.t2;
extrinsics = geom.betas;

%---loop through alongshore locations
num = -25:-1:-450;
for i = 1:length(num)
    
    m = R.slope(pindex,tindex);

    R2.z(1) = R.runup05(pindex,tindex);
    R2.z(2) = R.runup(pindex,tindex);
    R2.z(3) = R.runup95(pindex,tindex);
    R2.x(1) = interp1(profile.z,profile.x,R2.z(1));
    R2.x(2) = interp1(profile.z,profile.x,R2.z(2));
    R2.x(3) = interp1(profile.z,profile.x,R2.z(3));
    R2.y = num(i)*ones(size(R2.x)); %repmat(TWL.y,size(TWL.x,2),1);
    
    TWL.z(1) = R.twl05(pindex,tindex);
    TWL.z(2) = R.twl(pindex,tindex);
    TWL.z(3) = R.twl95(pindex,tindex);
    TWL.x(1) = interp1(profile.z,profile.x,TWL.z(1));
    TWL.x(2) = interp1(profile.z,profile.x,TWL.z(2));
    TWL.x(3) = interp1(profile.z,profile.x,TWL.z(3));
    TWL.y = num(i)*ones(size(TWL.x)); %repmat(TWL.y,size(TWL.x,2),1);
    
    %*** 4. world to image coords
    
    %---profile
    xyz = [profile.x,  profile.y, profile.z]; % may to transpose or not to get to work, also , or ;
    
    UV = round(xyz2DistUV(intrinsics, extrinsics, xyz));
    UV = reshape(UV,[],2);
    profile.u = UV(:,1);
    profile.v = UV(:,2);
    clear xyz UV
    
    %---R2
    xyz = [R2.x;  R2.y; R2.z]';
    UV = round(xyz2DistUV(intrinsics, extrinsics, xyz));
    UV = reshape(UV,[],2);
    R2.u(:,i) = UV(:,1);
    R2.v(:,i) = UV(:,2);
    clear xyz UV
    
    %---TWL
    xyz = [TWL.x;  TWL.y; TWL.z]';
    UV = round(xyz2DistUV(intrinsics, extrinsics, xyz));
    UV = reshape(UV,[],2);
    TWL.u(:,i) = UV(:,1);
    TWL.v(:,i) = UV(:,2);
%     clear xyz UV
end

figure(1);
hold on
% plot(profile.u,profile.v,'b-','linewidth',2);
% plot(TWL.u,TWL.v,'-r');
inan = find(isnan(TWL.u(1,:))==0);
PP = patch([TWL.u(1,end) TWL.u(3,end) TWL.u(3,inan(1)) TWL.u(1,inan(1))],...
    [TWL.v(1,end) TWL.v(3,end) TWL.v(3,inan(1)) TWL.v(1,inan(1))],'c');
set(PP,'FaceAlpha',.2); set(PP,'EdgeColor','none')
%hold on;  plot(TWL.u(2,:), TWL.v(2,:),'-b','linewidth',3);
hold on;  plot(TWL.u(2,:), TWL.v(2,:),'color',[15/255,90/255,185/255],'linewidth',3);
%hold on;  plot(R2.u(2,:), R2.v(2,:),'-c','linewidth',2);
% set(gca,'xticklabel','','yticklabel','');
ax = axis;
axis image; axis(ax);
%title([datestr(R.time(tindex)) 'GMT']);
print('-dpng', ['C:\Users\eswanson\OneDrive - DOI\Documents\GitHub\pixel2runup\forecasted',num2str(year), num2str(month), num2str(day), num2str(hr1), '.png'])

figure(2);
hold on
% plot(profile.u,profile.v,'b-','linewidth',2);
% plot(TWL.u,TWL.v,'-r');
inan = find(isnan(TWL.u(1,:))==0);
PP = patch([TWL.u(1,end) TWL.u(3,end) TWL.u(3,inan(1)) TWL.u(1,inan(1))],[TWL.v(1,end) TWL.v(3,end) TWL.v(3,inan(1)) TWL.v(1,inan(1))],'c');
set(PP,'FaceAlpha',.2); set(PP,'EdgeColor','none')
hold on;  plot(TWL.u(2,:), TWL.v(2,:),'color',[15/255,90/255,185/255],'linewidth',3);
%hold on;  plot(R2.u(2,:), R2.v(2,:),'-c','linewidth',3);
% set(gca,'xticklabel','','yticklabel','');
ax = axis;
axis image; axis(ax);
%title([datestr(R.time(tindex)) 'GMT']);
print('-dpng', 'C:\Users\eswanson\OneDrive - DOI\Documents\GitHub\pixel2runup\forecasted202291513_pt2.png')

% %*** 5. measured data
% load(['\\gs\StPetersburgFL-G\NACCH\Imagery\madbeach\runup\',num2str(year),'\',yearday,'_',monthName,'.',dayNum,'\',num2str(year),monthNum,dayNum,'_',num2str(hr1),'00_runup90_fromProfile.mat']);
% wl.z(1) = Runup.R2;
% wl.z(2) = Runup.TWL;
% wl.z(3) = Runup.tide;
% wl.x(1) = interp1(profile.z,profile.x,wl.z(1));
% wl.x(2) = interp1(profile.z,profile.x,wl.z(2));
% wl.x(3) = interp1(profile.z,profile.x,wl.z(3));
% 
% wl.x(4) = Runup.X2;
% wl.z(4) = interp1(profile.x,profile.z,wl.x(4));
% wl.x(5) = Runup.Xi2;
% wl.z(5) = interp1(profile.x,profile.z,wl.x(5));
% 
% wl.y = -90*ones(size(wl.x));
% xyz = [wl.x;  wl.y; wl.z]';
% UV = round(findUVnDOF(geom.betas,xyz,geom.meta.globals));
% UV = reshape(UV,[],2);
% wl.u = UV(:,1);
% wl.v = UV(:,2);
% 
% xyz_all = [profile.x, -90*ones(size(profile.x)), profile.z];
% % xyz_all = [profile.x, profile.y, profile.z];
% UV_all = round(findUVnDOF(geom.betas,xyz_all,geom.meta.globals));
% UV_all = reshape(UV_all,[],2);
% 
% figure(1)
% %plot(wl.u(1),wl.v(1),'bo','markersize',10,'markerfacecolor','b')
% %plot(UV_all(:,1),UV_all(:,2),'linewidth',1)
% plot(wl.u(2),wl.v(2),'go','markersize',10,'markerfacecolor','g')
% %plot(wl.u(3),wl.v(3),'mo','markersize',10,'markerfacecolor','m')
% 
% %plot(wl.u(4),wl.v(4),'bs','markersize',10);
% %plot(wl.u(5),wl.v(5),'bo','markersize',10);
% 
% 
% %legend('95% uncertainty band','Total Water Level Forecast','model R2','meas R2','meas TWL','meas X2','meas Xi2');%,'meas tide')
% %legend('','model TWL','model R2','meas R2','meas TWL','meas X2','meas Xi2','line 60');%,'meas tide')
% 
% %legend('95% uncertainty band','Total Water Level Forecast','Observed Total Water Level', 'runup transect')
% print('-dpng', ['C:\Users\mpalmsten\OneDrive - DOI\projects\TWL\MadeiraBeach\ForecastOnImage\SnapTWLForecastAndData',num2str(year), num2str(month), num2str(day), num2str(hr1), '.png'])

figure(2)
%plot(wl.u(1),wl.v(1),'bo','markersize',10,'markerfacecolor','b')
%plot(UV_all(:,1),UV_all(:,2),'linewidth',1, 'color','b')
plot(wl.u(2),wl.v(2),'go','markersize',10,'markerfacecolor','g')
%plot(wl.u(3),wl.v(3),'mo','markersize',10,'markerfacecolor','m')

%plot(wl.u(4),wl.v(4),'bs','markersize',10);
%plot(wl.u(5),wl.v(5),'bo','markersize',10);

foo = gca;
foo.XTick = '';
foo.YTick = '';
print('-dpng', ['C:\Users\mpalmsten\OneDrive - DOI\projects\TWL\MadeiraBeach\ForecastOnImage\TimexTWLForecastAndData',num2str(year), num2str(month), num2str(day), num2str(hr1), '.png'])

%legend('95% uncertainty band','total water level forecast','runup transect','observed total water level')
%print -dpng C:\Users\mpalmsten\Documents\presentations\StPeteScienceFest\2021\TWL2020091616.png
%legend('','model TWL','model R2','meas R2','meas TWL','meas X2','meas Xi2');%,'meas tide')

% keyboard
% if 1
%     saveas(figure(1),['\\gs\StPetersburgFL-G\NACCH\Imagery\madbeach\plots\snapshot_R2_',num2str(year),'_',monthNum,'_',dayNum,'_',num2str(hr1+1),'00'],'png')
% end

camera.x = geom.betas(1);
camera.y = geom.betas(2);
camera.z = geom.betas(3);

% [xiRunup,yiRunup,ziRunup]=computeRunupMapped2Topo(Runup.X2,-90,0,topo.x,topo.y,topo.z,camera.x,camera.y,camera.z);
[xiRunup,yiRunup,ziRunup]=computeRunupMapped2Topo(Runup.X2,-90,0,profile.x,profile.y,profile.z,camera.x,camera.y,camera.z);


figure;
subplot(2,1,1);
plot(profile.x,profile.z,'k','linewidth',2);
hold on
lh(1)=plot(Runup.X2,wl.z(4),'bs');
lh(2)=plot(Runup.Xi2,Runup.R2,'bo');
lh(3)=plot(wl.x(1),wl.z(1),'c*');
lh(4)=plot(wl.x(2),wl.z(2),'r*');
lh(5)=plot(R2.x(2), R2.z(2),'c^');
lh(6)=plot(TWL.x(2), TWL.z(2),'r^');
lh(7)=plot(interp1(profile.z,profile.x,Runup.param.R2),Runup.param.R2,'ms');
plot(xiRunup,ziRunup,'kp');
set(gca,'xdir','reverse');
legend(lh,'X2%','X2% (interp)','R2% (z interp to x-profile)','TWL (z interp to x-profile)',...
    'R2% (model)','TWL (model)','Stockdon param','location','southeast');
title([datestr(R.time(tindex)) 'GMT']);
xlabel('cross-shore distance');ylabel('elevation')

clear lh
subplot(2,1,2);
plot(profile.u,profile.v,'k','linewidth',2);
hold on
lh(1)=plot(wl.u(1),wl.v(1),'c*');
lh(2)=plot(wl.u(2),wl.v(2),'r*');
iloc=find(num==-90);
lh(3)=plot(R2.u(2,iloc), R2.v(2,iloc),'c^');
lh(4)=plot(TWL.u(2,iloc), TWL.v(2,iloc),'r^');
plot(wl.u(4),wl.v(4),'bs');
plot(wl.u(5),wl.v(5),'bo');
xlabel('x-position on image');ylabel('y-position on image')
legend(lh,'R2% (z interp to x-profile)','TWL (z interp to x-profile)','R2% (model)','TWL (model)','location','southeast');


