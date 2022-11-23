%Madeira Beach, FL - Plot all data together

%---ADCP (27.78897, -82.81229, )
%       datum: h_ADCP(NAVD88) = z_ADCP (bathy, NAVD88) + dp_ADCP (sens.info.cell1) + h_ADCP
%              tide (NAVD88) = h_ADCP (NAVD88) - namneam(h_ADCP (NAVD88))


%---NOAA station 8726724, at Clearwater Beach (NAVD88)

%---ADCIRC

%---TWL Model (MSL, convert to NAVD88)

clear,close all

%DIR = 'D:\CoastCams\madbeach';
DIR = '\\gs\StPetersburgFL-G\NACCH\Imagery\madbeach';

calc=0;
%--------------------------------------------------------------------------
if calc==1
%*** ADCP

%***altitude above bottom (of pressure transducer)
dp = [0.69, 0.68, 0.67, 0.66, 0.63, 0, 0.67, 0];

%***Daylight Savings Times offsets
ET.dnum   = [datenum(2017,3,12,2,0,0),  datenum(2017,11,5,2,0,0),   datenum(2018,3,11,2,0,0),   datenum(2018,11,4,2,0,0),   datenum(2019,3,10,2,0,0),   datenum(2019,11,3,2,0,0)];
ET.offset = [4,                         5,                          4,                          5,                          4,                          5];

DIR_adcp = '\\gs\stpetersburgfl-g\NACCH\Imagery\madbeach\waves\ADCP';

all.t=[];
all.t_local=[];
all.Hs=[];
all.HsSwell=[];
all.HsSea=[];
all.T01=[];
all.Tp=[];
all.Dm=[];
all.pd=[];
all.h=[];
all.temp=[];
all.deployment = [];


D = dironly([DIR_adcp '\20*_20*']);
for dd=1:length(D)
    if D(dd).isdir
        FILE = dironly([DIR_adcp '\' D(dd).name '\processed*\currents\*.mat']);
        if ~isempty(FILE)
            %---load data (processed/exported with TRDI Velocity software)
            [data]=plotData_TRDISentinelV([FILE.folder '\' FILE.name],[],0);
            
            %---convert local time to GMT, based on date/time of start of deployment
            ind = find(data.waves.dnum(1) >= ET.dnum,1,'last');
            
            all.t = [all.t; data.sens.dnum + datenum(0,0,0,ET.offset(ind),0,0)];
            all.t_local = [all.t_local; data.sens.dnum];
            all.Hs = [all.Hs; data.waves.Hs];
            all.HsSwell = [all.HsSwell; data.waves.HsSwell];
            all.HsSea = [all.HsSea; data.waves.HsSea];
            all.T01 = [all.T01; data.waves.T01];
            all.Tp = [all.Tp; data.waves.Tp];
            all.Dm = [all.Dm; data.waves.Dm];
            all.pd = [all.pd; data.sens.pd];
            %all.h = [all.h; data.sens.pd + data.info.cell1]; %h = p + dp (sensor height above bottom)
            all.h = [all.h; data.sens.pd + dp(dd)]; %h = p + dp (sensor height above bottom)
            all.temp = [all.temp; data.sens.t];
            
            all.deployment = [all.deployment; ones(size(data.sens.dnum)).*dd];
            all.GMToffset(dd) = ET.offset(ind);
            
            clear data
        end
        clear FILE
    end
end
ADCP_all=all;
clear all


%***compute tide in NAVD88 (using bathy data)
ADCP_all.tide = ADCP_all.h - nanmean(ADCP_all.h);

save([DIR_adcp '\allADCP_to20190913.mat'],'ADCP_all');

elseif calc==0
    load([DIR '\waves\ADCP\allADCP_to20190913.mat']);
    ADCP_all.tide(3024)=NaN;
end

% %***check h (=p + dp)
% figure;
% s1(1)=subplot(2,1,1);
% plot(ADCP_all.t,ADCP_all.h,'-b.');
% hold on
% plot(ADCP_all.t,ADCP_all.pd,'-r.');
% s1(2)=subplot(2,1,2);
% plot(ADCP_all.t,ADCP_all.h-ADCP_all.pd,'-k.');
% linkaxes(s1,'x');
% dynamicDateTicks(s1,'linked','mm/dd/yy');

%--------------------------------------------------------------------------
%---NOAA tide gauge (Clearwater Beach, NOAA Station 8726724)
%   --> times are GMT, tide data is NAVD88

data=load([DIR '\tides\NOAA_allTides.mat']);
NOAA.t = data.tide.t;
NOAA.tide = data.tide.measured;

NOAA.tidei = interp1(NOAA.t,NOAA.tide,ADCP_all.t);

% data=load('D:\CoastCams\madbeach\tides\allTides.mat');
% NOAA.t=data.tideAll.t(:);
% NOAA.tide=data.tideAll.z(:);
% 
% NOAA.tidei = interp1(NOAA.t,NOAA.tide,ADCP_all.t);


clear data
%--------------------------------------------------------------------------
%---ADCIRC tides (@ ADCP)

data=load([DIR '\tides\ADCIRC\Madeira_ADCP_tide_2017_2020.mat']); %time in GMT (assumed)
ADCIRC_ADCP.t=data.T;
ADCIRC_ADCP.tide=data.tid;

ADCIRC_ADCP.tidei = interp1(ADCIRC_ADCP.t,ADCIRC_ADCP.tide,ADCP_all.t);

clear data
%--------------------------------------------------------------------------
%---ADCIRC tides (@ 20-m)
data=load([DIR '\tides\ADCIRC\ADCIRC_all_Madeira_ADCPloc2.mat']); %time in GMT (assumed)
ADCIRC_20m.t=data.TT;
ADCIRC_20m.tide=data.zeta20;

ADCIRC_20m.tidei = interp1(ADCIRC_20m.t,ADCIRC_20m.tide,ADCP_all.t); %dates do not overlap

ADCIRC=data;

clear data
%--------------------------------------------------------------------------
%---TWL Model

MSL2NAVD88 = -0.096;

data=load('\\gs\stpetersburgfl-g\NACCH\Imagery\madbeach\runup\all_TWL_forecast.mat');

TWL.t=datenum({data.all_twl(:).forecastTime});
TWL.tide(:,1)=[[data.all_twl(:).twl]-[data.all_twl(:).runup]] + MSL2NAVD88; %converted to NAVD

[TWL.t,inds]=unique(TWL.t);
TWL.tide=TWL.tide(inds);
TWL.tidei=interp1(TWL.t,TWL.tide,ADCP_all.t);

TWL.Hs(:,1)=[data.all_twl(:).hs]; TWL.Hs=TWL.Hs(inds);
TWL.Tp(:,1)=[data.all_twl(:).tp]; TWL.Tp=TWL.Tp(inds);

clear data inds
%--------------------------------------------------------------------------
%***ADCP timeseries
figure;
ax(1)=subplot(5,1,1);
plot(ADCP_all.t,ADCP_all.HsSwell,'b');
hold on
plot(ADCP_all.t,ADCP_all.HsSea,'r');
plot(ADCP_all.t,ADCP_all.Hs,'k.','linewidth',2);
legend('swell','sea');
ylabel('H_s (m)');
set(gca,'xticklabel','');
xlim([ADCP_all.t(1) ADCP_all.t(end)]);

ax(2)=subplot(5,1,2);
plot(ADCP_all.t,ADCP_all.T01,'r','linewidth',2);
hold on
plot(ADCP_all.t,ADCP_all.Tp,'color',[0.5 0.5 0.5],'linewidth',2);
%plot(waves.dnum,waves.TpSwell,'b--','linewidth',2);
%plot(waves.dnum,waves.TpSea,'r--','linewidth',2);
legend('mean','peak');%,'peak-swell','peak-sea');
ylim([0 15]);
ylabel('T_m (s)');
set(gca,'xticklabel','');
xlim([ADCP_all.t(1) ADCP_all.t(end)]);

ax(3)=subplot(5,1,3);
plot(ADCP_all.t,ADCP_all.Dm,'g','linewidth',2);
hold on
% plot([waves.dnum(1) waves.dnum(end)],[40 40],'k--');
% plot([waves.dnum(1) waves.dnum(end)],[130 130],'k--');
% plot([waves.dnum(1) waves.dnum(end)],[220 220],'k--');
% plot([waves.dnum(1) waves.dnum(end)],[310 310],'k--');
%plot(waves.dnum,waves.Dp,'k--','linewidth',2);
%legend('mean','peak','peak-swell','peak-sea');
ylabel('\theta_m');
set(gca,'xticklabel','');
xlim([ADCP_all.t(1) ADCP_all.t(end)]);

ax(4)=subplot(5,1,4);
plot(ADCP_all.t,ADCP_all.h,'b','linewidth',2);
hold on
set(gca,'xticklabel','');
ylabel('h (m)');
xlim([ADCP_all.t(1) ADCP_all.t(end)]);

ax(5)=subplot(5,1,5);
plot(ADCP_all.t,ADCP_all.temp,'k','linewidth',2);
hold on
ylabel('temperature');
%datetick('x','mm/dd','keepticks')
xlabel('time');
linkaxes(ax,'x');
% dynamicDateTicks(ax,'linked','mm/dd/yy');
datetick('x')



%***model timeseries comparisons (no overlap with ADCIRC (nowcast) and ADCP data)
t1=find(TWL.t>=datenum(2020,1,1));
ta=find(ADCIRC_ADCP.t>=datenum(2020,1,1,0,0,0));

figure;
xx(1)=subplot(5,1,1);
%plot(ADCP_all.t,ADCP_all.Hs,'k.');
%hold on
plot(TWL.t(t1),TWL.Hs(t1),'-m.');
hold on
plot(ADCIRC.TT,ADCIRC.Hs,'-b.');
hold on
plot(ADCIRC.TT,ADCIRC.Hs20,'-r.');
%legend('ADCP','ADCIRC (5m)','ADCIRC (20m)');
legend('TWL (20m)','ADCIRC (5m, nowcast)','ADCIRC (20m, nowcast)');
title('model comparisons');
ylabel('H_s (m)'); 
set(gca,'xticklabel','');
%xlim([ADCP_all.t(1) ADCP_all.t(end)]);

xx(2)=subplot(5,1,2);
%plot(ADCP_all.t,ADCP_all.Tp,'-k.');
%hold on
plot(TWL.t(t1),TWL.Tp(t1),'-m.');
hold on
plot(ADCIRC.TT,ADCIRC.Tp,'-b.');
hold on
plot(ADCIRC.TT,ADCIRC.Tp20,'-r.');
ylim([0 15]);
ylabel('T_p (s)');
set(gca,'xticklabel','');
%xlim([ADCP_all.t(1) ADCP_all.t(end)]);

xx(3)=subplot(5,1,3);
%plot(ADCP_all.t,ADCP_all.Dm,'-k.','linewidth',2);
%hold on
plot(ADCIRC.TT,ADCIRC.Wdir,'-b.');
hold on
plot(ADCIRC.TT,ADCIRC.Wdir20,'-r.');
ylabel('\theta_m');
set(gca,'xticklabel','');
%xlim([ADCP_all.t(1) ADCP_all.t(end)]);

xx(4)=subplot(5,1,4);
plot(TWL.t(t1),TWL.tide(t1),'-m.');
hold on
plot(ADCIRC_ADCP.t(ta),ADCIRC_ADCP.tide(ta),'-k.'); %predited
hold on
plot(ADCIRC.TT,ADCIRC.zeta,'-b.'); %nowcast
plot(ADCIRC_20m.t,ADCIRC_20m.tide,'--r.'); %nowcast
legend('20m (TWL)','5m (predicted)','5m (nowcast)','20m (nowcast)');
ylabel('tide');
%datetick('x','mm/dd','keepticks')
%xlim([ADCP_all.t(1) ADCP_all.t(end)]);
linkaxes(xx,'x');


tidei_5m =interp1(ADCIRC_ADCP.t(ta),ADCIRC_ADCP.tide(ta),ADCIRC.TT);
tidei_TWL=interp1(TWL.t,TWL.tide,ADCIRC.TT);

xx(5)=subplot(5,1,5);
plot(ADCIRC.TT,ADCIRC.zeta-tidei_5m,'-k.'); %5m (nowcast) - 5m (predicted)
hold on
plot(ADCIRC.TT,ADCIRC.zeta20-tidei_5m,'--r.');       %20m (nowcast) - 5m (predicted)
plot(ADCIRC.TT,ADCIRC.zeta20-ADCIRC.zeta,'-b.'); %20m (nowcast) - 5m (nowcast)
plot(ADCIRC.TT,ADCIRC_20m.tide-tidei_TWL,'-m.');
legend('5m (nowcast) - 5m (predicted)','20m (nowcast) - 5m (predicted)','20m (nowcast) - 5m (nowcast)','20m (nowcast) - 20m (TWL)');
xlabel('time');
linkaxes(xx,'x');
% dynamicDateTicks(xx,'linked','mm/dd/yy');
datetick('x')



%---tides
figure;
axx(1)=subplot(3,1,1);
lh(2)=plot(ADCIRC_ADCP.t,ADCIRC_ADCP.tide,'-b.'); %predicted
hold on
%plot(ADCP_all.t,ADCIRC_ADCP.tidei,'--k.');
lh(3)=plot(ADCIRC_20m.t,ADCIRC_20m.tide,'-r.'); %nowcast
lh(4)=plot(TWL.t,TWL.tide,'-m.');
lh(5)=plot(NOAA.t,NOAA.tide,'-g.');
lh(1)=plot(ADCP_all.t,ADCP_all.tide,'-k.');
legend(lh,'ADCP (5m)','ADCIRC (5m, predicted)','ADCIRC (20m, nowcast)','TWL Model (20m)','NOAA');

axx(2)=subplot(3,1,2);
ll(2)=plot(ADCIRC_ADCP.t,ADCIRC_ADCP.tide,'-b.');
hold on
ll(1)=plot(ADCP_all.t,ADCP_all.tide,'-k.');
legend(ll,'ADCP (5m)','ADCIRC (5m, predicted)');

axx(3)=subplot(3,1,3);
lz(2)=plot(ADCP_all.t,ADCP_all.tide-TWL.tidei,'-m.');
hold on
lz(1)=plot(ADCP_all.t,ADCP_all.tide-ADCIRC_ADCP.tidei,'-k.');
legend(lz,'measured (ADCP, 5m) - predicted (ADCIRC, 5m)','measured (ADCP, 5m) - predicted (TWL, 20m)');
linkaxes(axx,'x');
% dynamicDateTicks(axx,'linked','mm/dd/yy');
datetick('x')


%---tides - linreg
[m_NOAA,b_NOAA,rsq_NOAA]=linreg(ADCP_all.tide,NOAA.tidei);
[m_5m,b_5m,rsq_5m]=linreg(ADCP_all.tide,ADCIRC_ADCP.tidei); %ADCIRC predicted

TWL.tidei_ADCIRC = interp1(TWL.t,TWL.tide,ADCIRC.TT);
[m_20m,b_20m,rsq_20m]=linreg(TWL.tidei_ADCIRC,ADCIRC_20m.tide);

figure;
subplot(2,2,1);
plot(ADCP_all.tide,ADCIRC_ADCP.tidei,'k.');
hold on
lg=plot(ADCP_all.tide,m_5m.*ADCP_all.tide+b_5m,'--k');
legend(lg,['m=' num2str(m_5m) ', r^2=' num2str(rsq_5m)]);
xlabel('ADCP (5m)');
ylabel('ADCIRC (5m, predicted)');
axis equal
axis square

subplot(2,2,2);
plot(ADCP_all.tide,NOAA.tidei,'k.');
hold on
lg=plot(ADCP_all.tide,m_NOAA.*ADCP_all.tide+b_NOAA,'--k');
legend(lg,['m=' num2str(m_NOAA) ', r^2=' num2str(rsq_NOAA)]);
xlabel('ADCP (5m)');
ylabel('NOAA (Clearwater pier)');
axis equal
axis square

subplot(2,2,3);
plot(TWL.tidei_ADCIRC,ADCIRC_20m.tide,'k.');
hold on
lg=plot(TWL.tidei_ADCIRC,m_20m.*TWL.tidei_ADCIRC+b_20m,'--k');
legend(lg,['m=' num2str(m_20m) ', r^2=' num2str(rsq_20m)]);
xlabel('TWL (20m)');
ylabel('ADCIRC (20m, nowcast)');
axis equal
axis square

tidei_ADCIRCpred=interp1(ADCIRC_ADCP.t,ADCIRC_ADCP.tide,ADCIRC.TT);
[mA5,bA5,rA5]=linreg(ADCIRC.zeta,tidei_ADCIRCpred);
[mA20,bA20,rA20]=linreg(ADCIRC.zeta,ADCIRC.zeta20);

subplot(2,2,4);
plot(ADCIRC.zeta,tidei_ADCIRCpred,'k.');
hold on
plot(ADCIRC.zeta,ADCIRC.zeta20,'b.');
lg(1)=plot(ADCIRC.zeta,mA5.*ADCIRC.zeta+bA5,'--k');
lg(2)=plot(ADCIRC.zeta,mA20.*ADCIRC.zeta+bA20,'--b');
legend(lg,['m=' num2str(mA5) ', r^2=' num2str(rA5)],['m=' num2str(mA20) ', r^2=' num2str(rA20)]);
xlabel('ADCIRC (5m, nowcast)');
ylabel({'ADCIRC (5m, predicted) - black';'ADCIRC (20m, nowcast) - blue'});
axis equal
axis square

return
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%---check EDT/EST to GMT conversions
crs=['r';'b';'g';'m';'c';'w';'y'];

figure;
subplot(2,1,1);
plot(ADCP_all.t_local,ADCP_all.tide,'-k.');
hold on
for tt=1:length(ET.dnum)
    plot([ET.dnum(tt) ET.dnum(tt)],[-2 2],'k--');
end
plot(ADCP_all.t_local(ADCP_all.deployment==1),ADCP_all.tide(ADCP_all.deployment==1),'ko','color',crs(1,:));
plot(ADCP_all.t_local(ADCP_all.deployment==2),ADCP_all.tide(ADCP_all.deployment==2),'ko','color',crs(2,:));
plot(ADCP_all.t_local(ADCP_all.deployment==3),ADCP_all.tide(ADCP_all.deployment==3),'ko','color',crs(3,:));
plot(ADCP_all.t_local(ADCP_all.deployment==4),ADCP_all.tide(ADCP_all.deployment==4),'ko','color',crs(4,:));
plot(ADCP_all.t_local(ADCP_all.deployment==5),ADCP_all.tide(ADCP_all.deployment==5),'ko','color',crs(5,:));
plot(ADCP_all.t_local(ADCP_all.deployment==6),ADCP_all.tide(ADCP_all.deployment==6),'ko','color',crs(6,:));
plot(ADCP_all.t_local(ADCP_all.deployment==7),ADCP_all.tide(ADCP_all.deployment==7),'ko','color',crs(7,:));
datetick('x','mm/dd/yy','keepticks');
xlabel('time (local)');


subplot(2,1,2);
plot(ADCIRC_ADCP.t,ADCIRC_ADCP.tide,'--k.','color',[0.6 0.6 0.6]);
hold on
plot(ADCP_all.t,ADCP_all.tide,'-k.');
hold on
%for tt=1:length(ET.dnum)
%    plot([ET.dnum(tt) ET.dnum(tt)],[-2 2],'k--');
%end
plot(ADCP_all.t(ADCP_all.deployment==1),ADCP_all.tide(ADCP_all.deployment==1),'ko','color',crs(1,:));
plot(ADCP_all.t(ADCP_all.deployment==2),ADCP_all.tide(ADCP_all.deployment==2),'ko','color',crs(2,:));
plot(ADCP_all.t(ADCP_all.deployment==3),ADCP_all.tide(ADCP_all.deployment==3),'ko','color',crs(3,:));
plot(ADCP_all.t(ADCP_all.deployment==4),ADCP_all.tide(ADCP_all.deployment==4),'ko','color',crs(4,:));
plot(ADCP_all.t(ADCP_all.deployment==5),ADCP_all.tide(ADCP_all.deployment==5),'ko','color',crs(5,:));
plot(ADCP_all.t(ADCP_all.deployment==6),ADCP_all.tide(ADCP_all.deployment==6),'ko','color',crs(6,:));
plot(ADCP_all.t(ADCP_all.deployment==7),ADCP_all.tide(ADCP_all.deployment==7),'ko','color',crs(7,:));

datetick('x','mm/dd/yy','keepticks');
xlabel('time (GMT)');
