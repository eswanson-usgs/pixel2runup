%% this file compares twl_output and the monthly mat files saved from the TWL viewer
% first need to run "load_TWL_forecasts.m" to get twl_output variable
% then run this to plot up Hs, etc. for the monthly file, then a plot to
% compare the two datasets.


clear R
load('\\gs\StPetersburgFL-G\NACCH\Projects\NWS\OutputArchive\WCOSS_test\monthFiles\3_2017runup.mat');

datevec(R.time(1))
datevec(R.time(end))

figure(1);clf
pcolor(R.Hs)
shading flat
colorbar
hline(200,'m');hline(330,'m');% for jan-march
% hline(420,'m');hline(900,'m');% for april-
% hline(1272,'k');hline(1535,'k');% for april-
xlabel('time')
ylabel('location')
title('Hs - April 2017')


% figure(2);clf
% plot(R.Lat_Update,'x')
% plot(R.longitude,'x')

cam.lat = 27.796216949206798;
cam.lon = -82.796102542950635;

pindex = find(min(abs(R.latitude-cam.lat))==abs(R.latitude-cam.lat))
% pindex = find(min(abs(R.Lat_Update-cam.lat))==abs(R.Lat_Update-cam.lat))

figure(3);clf
plot(R.time,R.Hs(pindex,:),'x-')
datetick('x',6)



%%
figure(4);clf
plot(date,vertcat(twl_output.twl),'o-')
hold on
plot(R.time,R.twl(pindex,:),'x--')
% plot(R.time,R.twl(pindex-1,:),'x--')
% plot(R.time,R.twl(pindex+1,:),'x--')
datetick('x',6,'keepticks')
xlabel('date')
ylabel('elevation')
title('TWL @ camera; monthly mat file vs. matlab query')
legend('matlab query','monthly mat file','monthly mat file pindex-1','monthly mat file pindex+1')