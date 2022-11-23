%% plot cross shore beach profile with estimated R2%
% also include tide elevation and TWL

%date to plot
year  = 2017;
month = 6;
day   = 1;

yearday=datenum(year,month,day)-datenum(year-1,12,31);
yearday=num2str(yearday,'%03d');
monthName=datestr([year,month,day,0,0,0],'mmm');
dayStr=num2str(day,'%2.2d');
dirName = [yearday '_' monthName '.' dayStr];

cd(['\\gs\StPetersburgFL-G\NACCH\Imagery\madbeach\runup\2017\',dirName])

runupfiles = dir([num2str(year),num2str(month,'%02g'),num2str(day,'%02g'),'*.mat'])

load(runupfiles(1).name)
figure(2);clf
plot(Runup.profile.x,Runup.profile.z,'b-','linewidth',3)
hold on

for ii=1:length(runupfiles)
    load(runupfiles(ii).name)
    
    tmp = Runup.R2-Runup.profile.z;
    [ind1,ind1]=min(abs(tmp));tmp(ind1) = NaN;
    [ind2,ind2]=min(abs(tmp));tmp1 = linspace(Runup.profile.z(ind1),Runup.profile.z(ind2),500);
    [ind3,ind3]=min(abs(Runup.R2-tmp1));tmp2 = linspace(Runup.profile.x(ind1),Runup.profile.x(ind2),500);
    R2x  = tmp2(ind3); clear tmp* ind*
    
    tmp = Runup.TWL-Runup.profile.z;
    [ind1,ind1]=min(abs(tmp));tmp(ind1) = NaN;
    [ind2,ind2]=min(abs(tmp));tmp1 = linspace(Runup.profile.z(ind1),Runup.profile.z(ind2),500);
    [ind3,ind3]=min(abs(Runup.TWL-tmp1));tmp2 = linspace(Runup.profile.x(ind1),Runup.profile.x(ind2),500);
    TWLx  = tmp2(ind3); clear tmp* ind*
    
    tmp = Runup.tide-Runup.profile.z;
    [ind1,ind1]=min(abs(tmp));tmp(ind1) = NaN;
    [ind2,ind2]=min(abs(tmp));tmp1 = linspace(Runup.profile.z(ind1),Runup.profile.z(ind2),500);
    [ind3,ind3]=min(abs(Runup.tide-tmp1));tmp2 = linspace(Runup.profile.x(ind1),Runup.profile.x(ind2),500);
    tidex  = tmp2(ind3); clear tmp* ind*
    
    figure(2)
    plot(R2x,Runup.R2,    'ko','markersize',10,'markerfacecolor','k')
    plot(TWLx,Runup.TWL,  'go','markersize',10,'markerfacecolor','g')
    plot(tidex,Runup.tide,'ro','markersize',10,'markerfacecolor','r')
    legend('Beach Profile','R2%','Total Water Level','Tide')
    title([runupfiles(ii).name(1:8),' - ',runupfiles(ii).name(10:13)])
    xlabel('cross-shore distance (m)');ylabel('elevation (m)');
    
    pause%(1)
end
plot(Runup.profile.x,0*Runup.profile.z,'k--','linewidth',1)
axis tight

