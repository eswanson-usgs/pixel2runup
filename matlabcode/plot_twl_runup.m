%% plot twl model vs runup timestack
clear 
clc
%% now plot comparisons of measured vs. forecast
cd \\gs\StPetersburgFL-G\NACCH\Imagery\madbeach

% run "load_all_TWL_forecasts.m" and "load_all_runup.m" to update .mat files

load('\\gs\StPetersburgFL-G\NACCH\Imagery\madbeach\runup\all_TWL_forecast.mat')
load('\\gs\StPetersburgFL-G\NACCH\Imagery\madbeach\runup\all_runup.mat')

    % remove runup files based off OLD elevation profile
 % this will always give an error ("index exceeds matrix dimensions b/c
    % we remove data...
    for ii=1:length(all_runup)-1
        if all_runup(ii).t==all_runup(ii+1).t
            all_runup(ii)=[];
            %will get an 'Index exceeds matrix dimensions' error at the
            %end, but its ok, continue to next section
        end
    end
    
  %%  
% !!!!! TWL model in MSL, measured data in NAVD88 -> subtract 0.087m!!!!!!
twl_output_twl        = vertcat(all_twl.twl)        -.087;
twl_output_twl05      = vertcat(all_twl.twl05)      -.087;
twl_output_twl95      = vertcat(all_twl.twl95)      -.087;
twl_output_runup      = vertcat(all_twl.runup)      ;%-.087;
twl_output_runup05    = vertcat(all_twl.runup05)    ;%-.087;
twl_output_runup95    = vertcat(all_twl.runup95)    ;%-.087;
% twl_output_tide    = vertcat(all_twl.tide)-.087;
twl_output_setup      = vertcat(all_twl.setup)      ;%-.087;
twl_output_swash      = vertcat(all_twl.swash)      ;%-.087;
twl_output_incSwash   = vertcat(all_twl.incSwash)   ;%-.087;
twl_output_infragSwash= vertcat(all_twl.infragSwash);%-.087;
% !!!!!!*******!!!!!!!******

% frcst_time = date';
frcst_time = datenum(vertcat(all_twl.forecastTime));
tmp = vertcat(all_runup.t);
tmp = datevec(tmp);
meas_time = datenum(tmp(:,1),tmp(:,2),tmp(:,3),tmp(:,4),0,0);
clear tmp
  dateStart = datenum(2017,01,21); 
  dateEnd   = datenum(2020,06,11);
md1 = find(meas_time>=dateStart,1,'first');
md2 = find(meas_time<=dateEnd,1,'last');
fd1 = find(frcst_time>=dateStart,1,'first');
fd2 = find(frcst_time<=dateEnd,1,'last');

%find dates/times where measured/forecasted are both available
lia1 = ismember(meas_time,frcst_time); % find frcst_time in meas_time
lia2 = ismember(frcst_time,meas_time); % find meas_time in frcst time
meas_d1=[]; frc_d1=[];
for ii = 1:length(frcst_time)
    if lia2(ii)==1
        frc_d1 = [frc_d1;ii];
    end
end
for ii = 1:length(meas_time)
    if lia1(ii)==1
        meas_d1 = [meas_d1;ii];
    end
end

if length(meas_d1) ~= length(frc_d1)
    error('corresponding dates dont match meas/frc')

end % will need to redo previous steps


%% plot TWL timeseries on top of each other

if(0) 
    figure(3);clf
    % plot(frcst_time(fd1:fd2),vertcat(twl_output(fd1:fd2).twl),'o-')
    plot(frcst_time(fd1:fd2),twl_output_twl(fd1:fd2),'o-','color',[.85 0.33 0.1],'markersize',10)%,'markerfacecolor',[.85 0.33 0.1])
    hold on
    plot(meas_time(md1:md2),vertcat(all_runup(md1:md2).TWL),'x','color',[0 0.45 0.74],'markersize',10)%,'markerfacecolor',[0 0.45 0.74])
    h=fill([frcst_time(fd1:fd2);flipud(frcst_time(fd1:fd2))],[twl_output_twl95(fd1:fd2);flipud(twl_output_twl05(fd1:fd2))],...
        [1 0.6 0.6],'edgecolor',[1 0.6 0.6]);
    uistack(h,'bottom')
    xlim([datenum(2017,6,16) datenum(2017,6,23)]) % define specific time range
    datetick('x',2,'keeplimits');
    title('Measured vs. Forecasted TWL')
    xlabel('Date');ylabel('Elevation (m, NAVD88)')
    legend('Forecast Uncertainty','Forecast','Measured')
    grid on
    
    %plot runup timeseries on top of each other
    figure(4);clf
    % plot(frcst_time(fd1:fd2),vertcat(twl_output(fd1:fd2).runup),'o-')
    plot(frcst_time(fd1:fd2),twl_output_runup(fd1:fd2),'o-','color',[.85 0.33 0.1],'markersize',10)
    hold on
    plot(meas_time(md1:md2),vertcat(all_runup(md1:md2).R2),'x','color',[0 0.45 0.74],'markersize',10)
    h=fill([frcst_time(fd1:fd2);flipud(frcst_time(fd1:fd2))],[twl_output_runup95(fd1:fd2);flipud(twl_output_runup05(fd1:fd2))],...
        [1 0.6 0.6],'edgecolor',[1 0.6 0.6]);
    uistack(h,'bottom')
    xlim([736862 736869])
    datetick('x',6,'keeplimits');
    title('Measured vs. Forecasted Runup')
    xlabel('Date');ylabel('Elevation (m, NAVD88)')
    legend('Forecast Uncertainty','Forecast','Measured')

    %plot tide timeseries on top of each other
    figure(5);clf
    % plot(frcst_time(fd1:fd2),vertcat(twl_output(fd1:fd2).tide),'o-')
    plot(frcst_time(fd1:fd2),twl_output_tide(fd1:fd2),'o-','color',[.85 0.33 0.1],'markersize',10)
    hold on
    plot(meas_time(md1:md2),vertcat(all_runup(md1:md2).tide),'x','color',[0 0.45 0.74],'markersize',10)
    % h=fill([frcst_time(fd1:fd2);flipud(frcst_time(fd1:fd2))],[twl_output_runup95(fd1:fd2);flipud(twl_output_runup05(fd1:fd2))],...
    %     [1 0.6 0.6],'edgecolor',[1 0.6 0.6]);
    % uistack(h,'bottom')
    xlim([736862 736869])
    datetick('x',6,'keeplimits');
    title('Measured vs. Forecasted Tide')
    xlabel('Date');ylabel('Elevation (m, NAVD88)')
    legend('Forecasted','Measured')
end

%% plot scatter of measured vs forecasted
%
% rrr = vertcat(twl_output(c).runup);
d1=datestr(all_runup(meas_d1(1)).t,23);
d2=datestr(all_runup(meas_d1(end)).t,23);

figure(6);clf
subplot(2,2,1) %runup
% % plot(vertcat(twl_output(frc_d1).twl),vertcat(all_runup(meas_d1).TWL),'x')
% % plot(twl_output_twl(frc_d1),vertcat(all_runup(meas_d1).TWL),'x')
% errorbar(vertcat(all_runup(meas_d1).TWL),twl_output_twl(frc_d1),...
%     abs(twl_output_twl(frc_d1)-twl_output_twl05(frc_d1)),abs(twl_output_twl(frc_d1)-twl_output_twl95(frc_d1)),...
%     '.','markersize',20)
% hold on
% plot(-0.6:.1:1.4,-0.6:.1:1.4,'k-')
% ylabel('Forecasted (m, NAVD88)');%xlabel('Measured (m, NAVD88)')
% title(['TWL - Dates: ',d1,' - ',d2])
% plot(vertcat(twl_output(frc_d1).runup),vertcat(all_runup(meas_d1).R2),'x')
% plot(vertcat(all_runup(meas_d1).R2),twl_output_runup(frc_d1),'.','markersize',20)
errorbar(vertcat(all_runup(meas_d1).R2),twl_output_runup(frc_d1),...
    abs(twl_output_runup(frc_d1)-twl_output_runup05(frc_d1)),abs(twl_output_runup(frc_d1)-twl_output_runup95(frc_d1)),...
    '.','markersize',20)
hold on
ax_neg=-0.5; ax_pos=1.5;
axis([ax_neg ax_pos ax_neg ax_pos])
axis square
plot(ax_neg:.1:ax_pos,ax_neg:.1:ax_pos,'k-')
ylabel('Forecasted (m, NAVD88)');%xlabel('Measured (m, NAVD88)')
title(['Runup Elevation - Dates: ',d1,' - ',d2])
[yest,slope,rsq,ypt,delta,rsq_sig] = quick_regress(vertcat(all_runup(meas_d1).R2),twl_output_runup(frc_d1));
text(-.35,1.25,['r^2 = ',num2str(round(rsq*100)/100)])

subplot(2,2,2) % setup
% plot(vertcat(twl_output(frc_d1).setup),vertcat(all_runup(meas_d1).setup),'x')
plot(vertcat(all_runup(meas_d1).setup),twl_output_setup(frc_d1),'.','markersize',20)
% errorbar(vertcat(all_runup(meas_d1).R2),twl_output_runup(frc_d1),'.','markersize',20)
hold on
ax_neg=-0.5; ax_pos=1;
axis([ax_neg ax_pos ax_neg ax_pos])
axis square
plot(ax_neg:.1:ax_pos,ax_neg:.1:ax_pos,'k-')
%ylabel('Forecasted (m, NAVD88)');%xlabel('Measured (m, NAVD88)')
title(['Setup Elevation'])
% do the following because there are some NaNs in all_runup.setup
tmp1=~isnan(vertcat(all_runup(meas_d1).setup));
tmp2=vertcat(all_runup(meas_d1(tmp1)).setup);
tmp3=twl_output_setup(frc_d1(tmp1));
% [yest,slope,rsq,ypt,delta,rsq_sig] = quick_regress(vertcat(all_runup(meas_d1).setup),twl_output_setup(frc_d1));
[yest,slope,rsq,ypt,delta,rsq_sig] = quick_regress(tmp2,tmp3);
text(-.35,0.85,['r^2 = ',num2str(round(rsq*100)/100)])

subplot(2,2,3) % incident swash
% plot(vertcat(all_runup(meas_d1).Sinc),twl_output_incSwash(frc_d1),'x')
plot(vertcat(all_runup(meas_d1).Sinc),twl_output_incSwash(frc_d1),'.','markersize',20)
% errorbar(vertcat(all_runup(meas_d1).Sinc),twl_output_incSwash(frc_d1),'.','markersize',20)
hold on
ax_neg=-0.2; ax_pos=1;
axis([ax_neg ax_pos ax_neg ax_pos])
axis square
plot(ax_neg:.1:ax_pos,ax_neg:.1:ax_pos,'k-')
ylabel('Forecasted (m, NAVD88)');xlabel('Measured (m, NAVD88)')
title(['Incident Swash Elevation'])
[yest,slope,rsq,ypt,delta,rsq_sig] = quick_regress(vertcat(all_runup(meas_d1).Sinc),twl_output_incSwash(frc_d1));
text(-0.1,0.9,['r^2 = ',num2str(round(rsq*100)/100)])

subplot(2,2,4) % infragravity swash
% plot(vertcat(all_runup(meas_d1).SIG),twl_output_infragSwash(frc_d1),'x')
plot(vertcat(all_runup(meas_d1).SIG),twl_output_infragSwash(frc_d1),'.','markersize',20)
% errorbar(vertcat(all_runup(meas_d1).SIG),twl_output_infragSwash(frc_d1),'.','markersize',20)
hold on
ax_neg=-0.25; ax_pos=1.5;
axis([ax_neg ax_pos ax_neg ax_pos])
axis square
plot(ax_neg:.1:ax_pos,ax_neg:.1:ax_pos,'k-')
xlabel('Measured (m, NAVD88)');%ylabel('Forecasted (m, NAVD88)');
title(['Infragravity Swash Elevation'])
[yest,slope,rsq,ypt,delta,rsq_sig] = quick_regress(vertcat(all_runup(meas_d1).SIG),twl_output_infragSwash(frc_d1));
text(-0.15,1.4,['r^2 = ',num2str(round(rsq*100)/100)])

% figure(8);clf
% % plot(vertcat(all_runup(meas_d1).tide),vertcat(twl_output(frc_d1).tide),'.','markersize',20)
% plot(vertcat(all_runup(meas_d1).tide),twl_output_tide(frc_d1),'.','markersize',20)
% hold on
% plot(-0.6:.1:0.8,-0.6:.1:0.8,'k-')
% ylabel('Forecasted');xlabel('Measured')
% title(['Measured vs. Forecast Tide Elevation (m, NAVD88)   -   Dates: ',d1,' - ',d2])


%% make timeseries plot of the variables

figure(9);clf
plot(vertcat(all_runup(meas_d1).t),vertcat(all_runup(meas_d1).R2),'.-','markersize',20)
hold on
plot(vertcat(all_runup(meas_d1).t),vertcat(all_runup(meas_d1).setup),'.-','markersize',20)
plot(vertcat(all_runup(meas_d1).t),vertcat(all_runup(meas_d1).Sinc),'.-','markersize',20)
plot(vertcat(all_runup(meas_d1).t),vertcat(all_runup(meas_d1).SIG),'.-','markersize',20)
hline(0,'k')
grid on

datetick('x','yyyymmdd')
xlabel('date');
ylabel('Elevation (m, NAVD88)');
legend('Runup','Setup','Incidenct Swash','Infragravity Swash')
title(['Measured Data - Dates: ',d1,' - ',d2])


figure(10);clf
errorbar(frcst_time(frc_d1),twl_output_runup(frc_d1),...
    abs(twl_output_runup(frc_d1)-twl_output_runup05(frc_d1)),abs(twl_output_runup(frc_d1)-twl_output_runup95(frc_d1)),...
    '.-','markersize',20)
hold on
plot(frcst_time(frc_d1),twl_output_setup(frc_d1),'.-','markersize',20)
plot(frcst_time(frc_d1),twl_output_incSwash(frc_d1),'.-','markersize',20)
plot(frcst_time(frc_d1),twl_output_infragSwash(frc_d1),'.-','markersize',20)
hline(0,'k')
grid on

datetick('x','yyyymmdd')
xlabel('date');
ylabel('Elevation (m, NAVD88)');
legend('Runup','Setup','Incidenct Swash','Infragravity Swash')
title(['Forecasted Data - Dates: ',d1,' - ',d2])

%% make timeseries of measured vs. modeled twl and r2
%make tmp variables to remove NaNs
% tmp1=~isnan(vertcat(all_runup(meas_d1).TWL));
% tmp2=vertcat(all_runup(meas_d1(tmp1)).TWL);
% tmp3=twl_output_twl(frc_d1(tmp1));
% bias is the mean of the difference between model and observations
bias_twl = nanmean(twl_output_twl(frc_d1)-vertcat(all_runup(meas_d1).TWL));
bias_twl_rel = nanmean(twl_output_twl(frc_d1)-vertcat(all_runup(meas_d1).TWL)/nanmean(twl_output_twl(frc_d1)));
bias_runup = nanmean(twl_output_runup(frc_d1)-vertcat(all_runup(meas_d1).R2));
bias_runup_rel = nanmean(twl_output_runup(frc_d1)-vertcat(all_runup(meas_d1).R2))/nanmean(twl_output_runup(frc_d1));
% scatter is Matlab 'rms' function on difference between observations and model
scatter_twl = rms(twl_output_twl(frc_d1)-vertcat(all_runup(meas_d1).TWL));
scatter_twl_index = rms(twl_output_twl(frc_d1)-vertcat(all_runup(meas_d1).TWL))/nanmean(vertcat(all_runup(meas_d1).TWL));
scatter_runup = rms(twl_output_runup(frc_d1)-vertcat(all_runup(meas_d1).R2));
scatter_runup_index = rms(twl_output_runup(frc_d1)-vertcat(all_runup(meas_d1).R2))/nanmean(vertcat(all_runup(meas_d1).R2));


figure(11);clf
plot(vertcat(all_runup(meas_d1).t),vertcat(all_runup(meas_d1).TWL),'-','markersize',20,'linewidth',2)
hold on
errorbar(frcst_time(frc_d1),twl_output_twl(frc_d1),...
    abs(twl_output_twl(frc_d1)-twl_output_twl05(frc_d1)),abs(twl_output_twl(frc_d1)-twl_output_twl95(frc_d1)),...
    '-','markersize',20,'linewidth',2)
h=hline(0,'k--');set(h,'linewidth',2);uistack(h,'bottom')
datetick('x','yyyy-mm-dd')
axis tight
xlim([all_runup(meas_d1(1)).t all_runup(meas_d1(end)).t])
legend('Measured','Modeled')
xlabel('Date (YYYY-MM-DD)')
ylabel('Elevation (m, NAVD88)')
title(['Madeira Beach Total Water Level (TWL): ',datestr(all_runup(meas_d1(1)).t),' - ',datestr(all_runup(meas_d1(end)).t)])

figure(12);clf
plot(vertcat(all_runup(meas_d1).t),vertcat(all_runup(meas_d1).TWL),'-','markersize',20,'linewidth',2)
hold on
plot(frcst_time(frc_d1),twl_output_twl(frc_d1),'-','markersize',20,'linewidth',2)
h=hline(0,'k--');set(h,'linewidth',2);uistack(h,'bottom')
datetick('x','yyyy-mm-dd')
axis tight
xlim([all_runup(meas_d1(1)).t all_runup(meas_d1(end)).t])
legend('Measured','Modeled')
xlabel('Date (YYYY-MM-DD)')
ylabel('Elevation (m, NAVD88)')
title(['Madeira Beach Total Water Level (TWL): ',datestr(all_runup(meas_d1(1)).t),' - ',datestr(all_runup(meas_d1(end)).t)])
text(datenum(2017,10,07)+3,1.4,'Hx Nate')
text(datenum(2018,05,27)+3,1.5,'TS Alberto')
text(datenum(2020,06,06)-80,1.1,'TS Cristobal')
text(datenum(2017,02,06),1.6,['model bias = ',num2str(round(bias_twl*100)/100)])
text(datenum(2017,02,06),1.525,['relative model bias = ',num2str(round(bias_twl_rel*100)/100)])
text(datenum(2017,02,06),1.45,['scatter = ',num2str(round(scatter_twl*100)/100)])
text(datenum(2017,02,06),1.375,['scatter index = ',num2str(round(scatter_twl_index*100)/100)])
text(datenum(2017,02,06),1.3,['n = ',num2str(length(vertcat(all_runup(meas_d1).TWL)))])

figure(13);clf
plot(vertcat(all_runup(meas_d1).t),vertcat(all_runup(meas_d1).R2),'-','markersize',20,'linewidth',2)
hold on
errorbar(frcst_time(frc_d1),twl_output_runup(frc_d1),...
    abs(twl_output_runup(frc_d1)-twl_output_runup05(frc_d1)),abs(twl_output_runup(frc_d1)-twl_output_runup95(frc_d1)),...
    '-','markersize',20,'linewidth',2)
h=hline(0,'k--');set(h,'linewidth',2);uistack(h,'bottom')
datetick('x','yyyy-mm-dd')
axis tight
xlim([all_runup(meas_d1(1)).t all_runup(meas_d1(end)).t])
legend('Measured','Modeled')
xlabel('Date (YYYY-MM-DD)')
ylabel('Elevation (m, NAVD88)')
title(['Madeira Beach Runup (R2%): ',datestr(all_runup(meas_d1(1)).t),' - ',datestr(all_runup(meas_d1(end)).t)])

figure(14);clf
plot(vertcat(all_runup(meas_d1).t),vertcat(all_runup(meas_d1).R2),'-','markersize',20,'linewidth',2)
hold on
plot(frcst_time(frc_d1),twl_output_runup(frc_d1),'-','markersize',20,'linewidth',2)
h=hline(0,'k--');set(h,'linewidth',2);uistack(h,'bottom')
datetick('x','yyyy-mm-dd')
axis tight
xlim([all_runup(meas_d1(1)).t all_runup(meas_d1(end)).t])
legend('Measured','Modeled')
xlabel('Date (YYYY-MM-DD)')
ylabel('Elevation (m, NAVD88)')
title(['Madeira Beach Runup (R2%): ',datestr(all_runup(meas_d1(1)).t),' - ',datestr(all_runup(meas_d1(end)).t)])
text(datenum(2017,02,06),1.3,['model bias = ',num2str(round(bias_runup*100)/100)])
text(datenum(2017,02,06),1.25,['relative model bias = ',num2str(round(bias_runup_rel*100)/100)])
text(datenum(2017,02,06),1.2,['scatter = ',num2str(round(scatter_runup*100)/100)])
text(datenum(2017,02,06),1.15,['scatter index = ',num2str(round(scatter_runup_index*100)/100)])
text(datenum(2017,02,06),1.1,['n = ',num2str(length(vertcat(all_runup(meas_d1).R2)))])


%% scatter of TWL and R2

d1=datestr(all_runup(meas_d1(1)).t,23);
d2=datestr(all_runup(meas_d1(end)).t,23);

figure(7);clf
subplot(2,2,1) %runup
% % plot(vertcat(twl_output(frc_d1).twl),vertcat(all_runup(meas_d1).TWL),'x')
% % plot(twl_output_twl(frc_d1),vertcat(all_runup(meas_d1).TWL),'x')
% errorbar(vertcat(all_runup(meas_d1).TWL),twl_output_twl(frc_d1),...
%     abs(twl_output_twl(frc_d1)-twl_output_twl05(frc_d1)),abs(twl_output_twl(frc_d1)-twl_output_twl95(frc_d1)),...
%     '.','markersize',20)
% hold on
% plot(-0.6:.1:1.4,-0.6:.1:1.4,'k-')
% ylabel('Forecasted (m, NAVD88)');%xlabel('Measured (m, NAVD88)')
% title(['TWL - Dates: ',d1,' - ',d2])
% plot(vertcat(twl_output(frc_d1).runup),vertcat(all_runup(meas_d1).R2),'x')
% plot(vertcat(all_runup(meas_d1).R2),twl_output_runup(frc_d1),'.','markersize',20)
errorbar(vertcat(all_runup(meas_d1).R2),twl_output_runup(frc_d1),...
    abs(twl_output_runup(frc_d1)-twl_output_runup05(frc_d1)),abs(twl_output_runup(frc_d1)-twl_output_runup95(frc_d1)),...
    '.','markersize',20)
hold on
ax_neg=-1; ax_pos=1.9;
axis([ax_neg ax_pos ax_neg ax_pos])
axis square
plot(ax_neg:.1:ax_pos,ax_neg:.1:ax_pos,'k-')
xlabel('Measured (m, NAVD88)');%ylabel('Forecasted (m, NAVD88)');
ylabel('Forecasted (m, NAVD88)');%xlabel('Measured (m, NAVD88)')
title(['Runup Elevation - Dates: ',d1,' - ',d2])
[yest,slope,rsq,ypt,delta,rsq_sig] = quick_regress(vertcat(all_runup(meas_d1).R2),twl_output_runup(frc_d1));
text(-.5,1.5,['r^2 = ',num2str(round(rsq*100)/100)])
display('Runup Elevation')
slope,rsq,ypt,delta,rsq_sig


subplot(2,2,2) % TWL
% plot(vertcat(twl_output(frc_d1).setup),vertcat(all_runup(meas_d1).setup),'x')
errorbar(vertcat(all_runup(meas_d1).TWL),twl_output_twl(frc_d1),...
    abs(twl_output_twl(frc_d1)-twl_output_twl05(frc_d1)),abs(twl_output_twl(frc_d1)-twl_output_twl95(frc_d1)),...
    '.','markersize',20)
% errorbar(vertcat(all_runup(meas_d1).R2),twl_output_runup(frc_d1),'.','markersize',20)
hold on
% ax_neg=-0.7; ax_pos=1.7;
axis([ax_neg ax_pos ax_neg ax_pos])
axis square
plot(ax_neg:.1:ax_pos,ax_neg:.1:ax_pos,'k-')
ylabel('Forecasted (m, NAVD88)');%xlabel('Measured (m, NAVD88)')
xlabel('Measured (m, NAVD88)');%ylabel('Forecasted (m, NAVD88)');
title(['TWL Elevation'])
tmp=isnan(vertcat(all_runup(meas_d1).TWL));
tmp1=find(tmp==1);
for ii=1:length(tmp1)
    tmp2=find(meas_d1==tmp1(ii));
    meas_d1(tmp2)=[];
    frc_d1(tmp1(ii)-ii+1)=[];
end
[yest,slope,rsq,ypt,delta,rsq_sig] = quick_regress(vertcat(all_runup(meas_d1).TWL),twl_output_twl(frc_d1));
text(-.5,1.5,['r^2 = ',num2str(round(rsq*100)/100)])
display('TWL Elevation')
slope,rsq,ypt,delta,rsq_sig

subplot(2,2,3) % TWL vs. tide/surge
plot(vertcat(all_runup(meas_d1).TWL),twl_output_twl(frc_d1)-twl_output_runup(frc_d1),'x')
% plot(vertcat(twl_output_tide(frc_d1)+twl_output_surge(frc_d1)),vertcat(all_runup(meas_d1).twl),'x')
% errorbar(vertcat(all_runup(meas_d1).TWL),twl_output_twl(frc_d1)-twl_output_runup(frc_d1),...
%     abs(twl_output_twl(frc_d1)-twl_output_runup(frc_d1)-twl_output_twl05(frc_d1)+twl_output_runup05(frc_d1)),...
%     abs(twl_output_twl(frc_d1)-twl_output_runup(frc_d1)-twl_output_twl95(frc_d1)+twl_output_runup95(frc_d1)),...
%     '.','markersize',20)
% errorbar(vertcat(all_runup(meas_d1).TWL),twl_output_twl(frc_d1),...
%     abs(twl_output_twl(frc_d1)-twl_output_twl05(frc_d1)),abs(twl_output_twl(frc_d1)-twl_output_twl95(frc_d1)),...
%     '.','markersize',20)
% errorbar(vertcat(all_runup(meas_d1).R2),twl_output_runup(frc_d1),'.','markersize',20)
hold on
% ax_neg=-0.7; ax_pos=1.7;
axis([ax_neg ax_pos ax_neg ax_pos])
axis square
plot(ax_neg:.1:ax_pos,ax_neg:.1:ax_pos,'k-')
ylabel('Forecasted Tide+Surge (m, NAVD88)');%xlabel('Measured (m, NAVD88)')
xlabel('Measured TWL (m, NAVD88)');%ylabel('Forecasted (m, NAVD88)');
title(['Measured TWL vs Modeled Tide+Surge'])
[yest,slope,rsq,ypt,delta,rsq_sig] = quick_regress(vertcat(all_runup(meas_d1).TWL),(twl_output_twl(frc_d1)-twl_output_runup(frc_d1)));
text(-.5,1.5,['r^2 = ',num2str(round(rsq*100)/100)])
display('Measured TWL vs Modeled Tide+Surge')
slope,rsq,ypt,delta,rsq_sig


%% scatter of TWL only

d1=datestr(all_runup(meas_d1(1)).t,23);
d2=datestr(all_runup(meas_d1(end)).t,23);

figure(8);clf
% TWL
errorbar(vertcat(all_runup(meas_d1).TWL),twl_output_twl(frc_d1),...
    abs(twl_output_twl(frc_d1)-twl_output_twl05(frc_d1)),abs(twl_output_twl(frc_d1)-twl_output_twl95(frc_d1)),...
    '.','markersize',20)
hold on
ax_neg=-1; ax_pos=1.9;
axis([ax_neg ax_pos ax_neg ax_pos])
plot(ax_neg:.1:ax_pos,ax_neg:.1:ax_pos,'k-')
ylabel('Forecasted (m, NAVD88)');%xlabel('Measured (m, NAVD88)')
xlabel('Measured (m, NAVD88)');%ylabel('Forecasted (m, NAVD88)');
title(['TWL Elevation - Dates: ',d1,' - ',d2])
axis square
tmp=isnan(vertcat(all_runup(meas_d1).TWL));
tmp1=find(tmp==1);
for ii=1:length(tmp1)
    tmp2=find(meas_d1==tmp1(ii));
    meas_d1(tmp2)=[];
    frc_d1(tmp1(ii)-ii+1)=[];
end
[yest,slope,rsq,ypt,delta,rsq_sig] = quick_regress(vertcat(all_runup(meas_d1).TWL),twl_output_twl(frc_d1),vertcat(all_runup(meas_d1).TWL));
text(-.5,1.5,['r^2 = ',num2str(round(rsq*100)/100)])
display('TWL Elevation')
slope,rsq,ypt,delta,rsq_sig

