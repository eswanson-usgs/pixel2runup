%Plot runup/TWL results

clear,close all

yd1 = date_to_yearday(2017,4,1);
yd2 = date_to_yearday(2017,4,18);
yds=[yd1:1:yd2];

t=[];
X2=[];
Xi2=[];
R2=[];
TWL=[];
tide=[];
R2param=[];
m=[];
Hs0=[];
Z2=[];
Zi2=[];
U2=[];
Ri2=[];
xRi2=[];

geom = load('\\gs\StPetersburgFL-G\NACCH\Imagery\madbeach\calibration\extrinsic\20170217\geomFile_20170217.mat');
camera.x = geom.betas(1);
camera.y = geom.betas(2);
camera.z = geom.betas(3);
load('\\gs\StPetersburgFL-G\NACCH\Imagery\madbeach\surveys\walking\20170217\20170217Local.mat','topo');

for ii=1:length(yds)
    [yr,mn,dy]=yearday_to_date(2017,yds(ii));
    
    R=dir(['\\gs\StPetersburgFL-G\NACCH\Imagery\madbeach\runup\2017\' num2str(yds(ii),'%3.3d') '_' monthstr(month(datenum(yr,mn,dy)),'mmm') '.' num2str(dy,'%2.2d') '\' num2str(yr) '*.mat']);
    
    for rr=1:length(R)
        load(['\\gs\StPetersburgFL-G\NACCH\Imagery\madbeach\runup\2017\' num2str(yds(ii),'%3.3d') '_' monthstr(month(datenum(yr,mn,dy)),'mmm') '.' num2str(dy,'%2.2d') '\' R(rr).name])
        
        t=[t; Runup.t];
        X2=[X2; Runup.X2];
        Xi2=[Xi2; Runup.Xi2];
        R2=[R2; Runup.R2];
        TWL=[TWL; Runup.TWL];
        tide=[tide; Runup.tide];
        R2param=[R2param; Runup.param.R2];
        m=[m; Runup.slope];
        Hs0=[Hs0; Runup.Hs0];
        good=find(isnan(Runup.profile.z)==0);
        Z2=[Z2; interp1(Runup.profile.x(good),Runup.profile.z(good),Runup.X2)];
        Zi2=[Zi2; interp1(Runup.profile.x(good),Runup.profile.z(good),Runup.Xi2)];
        U2=[U2; interp1(Runup.profile.z(good),Runup.profile.x(good),Runup.R2)];

        [xiRunup,yiRunup,ziRunup]=computeRunupMapped2Topo(Runup.X2,-90,0,topo.x,topo.y,topo.z,camera.x,camera.y,camera.z);
        Ri2=[Ri2; ziRunup];
        xRi2=[xRi2; xiRunup];
        
        clear Runup
    end
    clear R yr mn dy
end

figure;
plot(X2,Xi2,'.');
hold on
plot(X2,xRi2,'g.');
plot([36 50],[36 50],'k');
xlabel('X2%');
ylabel('X2% - interpolated to topo');

% figure;
% plot(R2,Hs0,'.');
% 
% figure;
% plot(R2,tide,'.');
% 
% figure;
% plot(R2,m,'.');

figure;
plot(R2,R2param,'.');
hold on
plot([0 1],[0 1],'k');
xlabel('R2% - measured');
ylabel('R2% - Stockdon');

figure;
plot(R2,Z2,'.');
hold on
plot(R2,Zi2,'r.');
plot(R2,Ri2,'g.')
plot([0 1],[0 1],'k');
xlabel('R2%');
ylabel('X2% interpolated to profile');

figure;
plot(X2,U2,'.');
hold on
plot(Xi2,U2,'r.');
plot(xRi2,U2,'g.');
plot([36 50],[36 50],'k');
xlabel('X2%');
ylabel('R2% interpolated to profile');

figure;
plot(t,tide,'ro');
hold on
for tt=1:18
    load(['\\gs\StPetersburgFL-G\NACCH\Imagery\madbeach\tides\tides_201704' num2str(tt,'%2.2d')]);
    
    plot(tides.t,tides.measured,'-k.');
    
    clear tides
end
datetick('x','mm/dd');
xlim([datenum(2017,4,1) datenum(2017,4,19)]);
    
    
    
    
    