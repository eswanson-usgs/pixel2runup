%Plot TWL estimates on raw camera image

clear,close all
addpath(genpath('C:\Users\jwlong\Documents\MATLAB\m-cmg\branches\jbrown\CIRN\CIL\UAV-Processing-Toolbox'))


year = 2017;
month = 3;  monthNum=num2str(month,'%02d');
day  = 31;  dayNum=num2str(day,'%02d');
hr  = 19;
yearday=datenum(year,month,day)-datenum(year-1,12,31);yearday=num2str(yearday,'%03d');
monthName=datestr([year,month,day,0,0,0],'mmm');
myepoch = num2str(datenum2epoch(datenum(year,month,day,hr,0,0)));
[~,dayName] = weekday(datenum(year,month,day));


%*** 1. image data
load(['\\gs\StPetersburgFL-G\NACCH\Archive\Data\2016\2016-363-DD_20161028\madbeach\2017\c1\',yearday,'_',monthName,'.',dayNum,'\',myepoch,'.',dayName,'.',monthName,'.',dayNum,'_',num2str(hr),'_00_00.GMT.',num2str(year),'.madbeach.c1.snap.mat']);
geom = load('\\gs\StPetersburgFL-G\NACCH\Imagery\madbeach\calibration\extrinsic\20170217\geomFile_20170217.mat');

figure(1);
imagesc(snap);

%*** 2. topo data
load('\\gs\StPetersburgFL-G\NACCH\Imagery\madbeach\surveys\walking\20170217\20170217Local.mat','topo');

yindex = find(min(abs(topo.y(:,1)-(-90)))==abs(topo.y(:,1)-(-90)));
profile.x = topo.x(yindex,:);
profile.z = topo.z(yindex,:);
good = find(isnan(profile.z)==0);
profile.x = profile.x(good);
profile.z = profile.z(good);
profile.y = -90.*ones(size(profile.x));

% figure(10);
% plot(profile.x,profile.z,'k');


%*** 3. model data
load(['\\gs\StPetersburgFL-G\NACCH\Projects\NWS\OutputArchive\WCOSS_test\monthFiles\',num2str(month),'_2017runup.mat']);

cam.lat = 27.796216949206798;
cam.lon = -82.796102542950635;

if month < 4 % the grid/latitude changed from March to April
pindex = find(min(abs(R.latitude-cam.lat))==abs(R.latitude-cam.lat));
elseif month >=4
pindex = find(min(abs(R.Lat_Update-cam.lat))==abs(R.Lat_Update-cam.lat));
end
tindex = find(R.time>datenum(year,month,day,hr,-30,0) & R.time<datenum(year,month,day,hr,30,0));

% added a loop to 
num = -25:-1:-450;
for i = 1:length(num)
m = R.slope(pindex,tindex);
TWL.z(1) = R.twl05(pindex,tindex);
TWL.z(2) = R.twl(pindex,tindex);
TWL.z(3) = R.twl95(pindex,tindex);
TWL.x(1) = interp1(profile.z,profile.x,TWL.z(1));
TWL.x(2) = interp1(profile.z,profile.x,TWL.z(2));
TWL.x(3) = interp1(profile.z,profile.x,TWL.z(3));
TWL.y = num(i)*ones(size(TWL.x)); %repmat(TWL.y,size(TWL.x,2),1);

% figure(10);
% hold on
% plot(TWL.x,TWL.z,'r*');

%*** 4. world to image coords

%---TWL
xyz = [TWL.x;  TWL.y; TWL.z]';
UV = round(findUVnDOF(geom.betas,xyz,geom.meta.globals));
UV = reshape(UV,[],2);
TWL.u(:,i) = UV(:,1);
TWL.v(:,i) = UV(:,2);
clear xyz UV

%---profile
xyz = [profile.x;  profile.y; profile.z]';
UV = round(findUVnDOF(geom.betas,xyz,geom.meta.globals));
UV = reshape(UV,[],2);
profile.u = UV(:,1);
profile.v = UV(:,2);
clear xyz UV

% now do runup
R2.z(1) = R.runup05(pindex,tindex);
R2.z(2) = R.runup(pindex,tindex);
R2.z(3) = R.runup95(pindex,tindex);
R2.x(1) = interp1(profile.z,profile.x,R2.z(1));
R2.x(2) = interp1(profile.z,profile.x,R2.z(2));
R2.x(3) = interp1(profile.z,profile.x,R2.z(3));
R2.y = num(i)*ones(size(R2.x)); %repmat(TWL.y,size(TWL.x,2),1);

% figure(10);
% hold on
% plot(TWL.x,TWL.z,'r*');

%*** 4. world to image coords

%---R2
xyz = [R2.x;  R2.y; R2.z]';
UV = round(findUVnDOF(geom.betas,xyz,geom.meta.globals));
UV = reshape(UV,[],2);
R2.u(:,i) = UV(:,1);
R2.v(:,i) = UV(:,2);
clear xyz UV

end

figure(1);
hold on
% plot(profile.u,profile.v,'b-','linewidth',2);
% plot(TWL.u,TWL.v,'-r');
inan = find(isnan(TWL.u(1,:))==0);
PP = patch([TWL.u(1,end) TWL.u(3,end) TWL.u(3,inan(1)) TWL.u(1,inan(1))],[TWL.v(1,end) TWL.v(3,end) TWL.v(3,inan(1)) TWL.v(1,inan(1))],'r');
set(PP,'FaceAlpha',.2); set(PP,'EdgeColor','none')
hold on;  plot(TWL.u(2,:), TWL.v(2,:),'-r','linewidth',3);  
hold on;  plot(R2.u(2,:), R2.v(2,:),'-c','linewidth',3);  
% set(gca,'xticklabel','','yticklabel','');
ax = axis; 
axis image; axis(ax);
title([datestr(R.time(tindex)) 'GMT']);

%now do measured runup from runup timestacks
    load(['\\gs\StPetersburgFL-G\NACCH\Imagery\madbeach\runup\',num2str(year),'\',yearday,'_',monthName,'.',dayNum,'\',num2str(year),monthNum,dayNum,'_',num2str(hr),'00_runup90.mat']);
    wl.z(1) = Runup.R2;
    wl.z(2) = Runup.TWL;
    wl.z(3) = Runup.tide;
    wl.x(1) = interp1(profile.z,profile.x,wl.z(1));
    wl.x(2) = interp1(profile.z,profile.x,wl.z(2));
    wl.x(3) = interp1(profile.z,profile.x,wl.z(3));
    wl.y = -90*ones(size(wl.x));;
    xyz = [wl.x;  wl.y; wl.z]';
    UV = round(findUVnDOF(geom.betas,xyz,geom.meta.globals));
    UV = reshape(UV,[],2);
    wl.u = UV(:,1);
    wl.v = UV(:,2);
    
    plot(wl.u(1),wl.v(1),'bo','markersize',10,'markerfacecolor','b')
    plot(wl.u(2),wl.v(2),'go','markersize',10,'markerfacecolor','g')
%     plot(wl.u(3),wl.v(3),'mo','markersize',10,'markerfacecolor','m')

legend('','model TWL','model R2','meas R2','meas TWL')%,'meas tide')

% keyboard
% if 1
%     saveas(figure(1),['\\gs\StPetersburgFL-G\NACCH\Imagery\madbeach\plots\snapshot_R2_',num2str(year),'_',monthNum,'_',dayNum,'_',num2str(hr),'00'],'png')
% end


figure;
subplot(2,1,1);
plot(profile.x,profile.z,'k','linewidth',2);
hold on
lh(1)=plot(Runup.X2,Runup.R2,'bs');
lh(2)=plot(Runup.Xi2,Runup.R2,'bo');
lh(3)=plot(wl.x(1),wl.z(1),'c*');
lh(4)=plot(wl.x(2),wl.z(2),'r*');
set(gca,'xdir','reverse');
legend(lh,'X2%','X2% (interp)','R2% (z interp to x-profile)','TWL (z interp to x-profile)');
title([datestr(R.time(tindex)) 'GMT']);

subplot(2,1,2);
plot(profile.u,profile.v,'k','linewidth',2);
hold on
lh(1)=plot(wl.u(1),wl.v(1),'c*');
lh(2)=plot(wl.u(2),wl.v(2),'r*');
iloc=find(num==-90);
lh(3)=plot(R2.u(2,iloc), R2.v(2,iloc),'c^');  
lh(4)=plot(TWL.u(2,iloc), TWL.v(2,iloc),'r^');  
legend(lh,'R2% (z interp to x-profile)','TWL (z interp to x-profile)','R2% (model)','TWL (model)');


