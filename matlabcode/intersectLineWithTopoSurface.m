%% Project runup line onto topo surface
% written by Meg Palmsten, USGS
% adapted by Justin Birchler, USGS for Madeira Beach, FL
cd '\\gs\stpetersburgfl-g\NACCH\Imagery\madbeach\'
camera_number = 'c2';
line_number = '52';
survey_date = '20220204';
addpath matlab
addpath surveys\walking

%  addpath C:\Users\mpalmsten\Documents\MATLAB\jbrown\RUNUP
addpath('C:\Users\jbirchler\OneDrive - DOI\Documents\stpete_branches\jbrown\RUNUP')

% load geometry etc
load(['\\gs\stpetersburgfl-g\NACCH\Imagery\madbeach\calibration\extrinsic\20220128\geomFile_',camera_number,'.mat'])
cams.x = betas(1);cams.y = betas(2);cams.z = betas(3);
% load('\\gs\stpetersburgfl-g\NACCH\Data\Molokai\geometry\geom\molokai.c2.r.mat')
% load('\\gs\stpetersburgfl-g\NACCH\Data\Molokai\geometry\geom\1529974978.c2.timex.molokai.meta.mat')
% load('\\gs\stpetersburgfl-g\NACCH\Data\Molokai\geometry\geom\moloGCPs_20180626_7lids_LOCAL.mat')
% load('\\gs\stpetersburgfl-g\NACCH\Data\Molokai\geometry\geom\moloGCPs_20180626_7lids_UTM.mat')
% 
% % load and set up topo
% topoFile.date = [2018 06 25];
% topoFile.lineno = [03];
% topoFile.path = ['\\gs\stpetersburgfl-g\NACCH\Data\Molokai\geometry\geom\'];
% [topoFile] = topoProfile(topoFile);

% load survey data
surveyData = dlmread(['.\surveys\walking\',survey_date,'\line',line_number,'.xyz']);

[Xsurvey,Ysurvey]=coordSys_madbeach(surveyData(:,1),surveyData(:,2));
Zsurvey = surveyData(:,3);
y_position = mean(Ysurvey);
xmin = -10; xmax = 100;
ymin = y_position-15; ymax = y_position+15;
% 
% % find PIX instance
% pixline.ind=find(r.cams.stacklineNum==3);
% pixline.u=r.cams.U(pixline.ind);
% pixline.v=r.cams.V(pixline.ind);
pixline.x=xmin:0.1:xmax;
pixline.y=ones(size(pixline.x))*y_position;
% pixline.z=r.cams.XYZ(pixline.ind,3);
pixline.z=ones(size(pixline.x))*0;
%
figure; 
hold on;grid on;box on;
plot(pixline.x,pixline.y,'--k')
scatter(Xsurvey,Ysurvey,50,Zsurvey,'filled')
set(gca,'ydir','reverse','xdir','reverse');colorbar
xlim([xmin xmax]); xlabel('cross-shore (x)')
ylim([ymin ymax]); ylabel('alongshore (y)')
title('surveyData, local coordinates')
caxis([-1.5 1.5])

% XYZ are GPS Survey data - use line data nearest runup line
X  = Xsurvey;
V  = Zsurvey;
Xq = (xmin:.1:xmax)';
Vq = interp1(X,V,Xq);
zF = Vq; % set zF to be output from interp1

[xSurfFull,ySurfFull] = meshgrid(xmin:.1:xmax,ymin:.1:ymax);
zFFull = ones(size(xSurfFull)).*zF';
zFFull(zFFull>max(Zsurvey))=max(Zsurvey);

%organize interpolated survey for use in processing
surveyInterp.x = xSurfFull(:);
surveyInterp.y = ySurfFull(:);
surveyInterp.z = zFFull(:);

figure; scatter(surveyInterp.x,surveyInterp.y,10,surveyInterp.z,'filled')
grid on;box on;set(gca,'ydir','reverse','xdir','reverse');colorbar
xlim([xmin xmax]); xlabel('cross-shore (x)')
ylim([ymin ymax]); ylabel('alongshore (y)')
title('surveyData, alongshore uniform, local coordinates')
caxis([-1.5 1.5])

%save interpolated beach topography
% save \\gs\stpetersburgfl-g\NACCH\Data\Molokai\geometry\topoSurvey.mat surveyInterp

% interpolate elevation of z=0 runup line to profile
% zTest = interp1(topoFile.profile.x,topoFile.profile.z,pixline.x);
zTest = interp1(Xsurvey,Zsurvey,pixline.x);
% zTest = Zsurvey;

% [xiRunup,yiRunup,ziRunup]=computeRunupMapped2Topo(pixline.x,pixline.y,pixline.z,xSurfFull,ySurfFull,zFFull,r.cams.x,r.cams.y,r.cams.z);
[xiRunup,yiRunup,ziRunup]=computeRunupMapped2Topo(pixline.x,pixline.y,pixline.z,xSurfFull,ySurfFull,zFFull,cams.x,cams.y,cams.z);
% [xiRunup,yiRunup,ziRunup]=computeRunupMapped2Topo(Xsurvey,Ysurvey,Zsurvey,xSurfFull,ySurfFull,zFFull,cams.x,cams.y,cams.z);
tmp1 = find(xiRunup>max(Xsurvey));
tmp2 = find(xiRunup<min(Xsurvey));
tmp4 = find(ziRunup<min(min(zFFull))); % had issue with large negative z values...
xiRunup(tmp1)=NaN;xiRunup(tmp2)=NaN;xiRunup(tmp4)=NaN;
yiRunup(tmp1)=NaN;yiRunup(tmp2)=NaN;yiRunup(tmp4)=NaN;
ziRunup(tmp1)=NaN;ziRunup(tmp2)=NaN;ziRunup(tmp4)=NaN;
[tmp3 ind3] = max(xiRunup); % find last non-Nan value

% plot first point in stack projected onto the beach
figure
hS = surf(xSurfFull,ySurfFull,zFFull,'faceColor','flat','EdgeColor','none')
hold on; grid on; box on
plot3(cams.x, cams.y, cams.z,'o')
text(cams.x+0.25, cams.y+0.25, cams.z+0.25,camera_number)
plot3(Xsurvey(1), Ysurvey(1), Zsurvey(1),'.b')
plot3([cams.x Xsurvey(1)], [cams.y Ysurvey(1)], [cams.z Zsurvey(1)],'-b')
plot3([cams.x xiRunup(ind3)], [cams.y yiRunup(ind3)], [cams.z ziRunup(ind3)],'-b')
plot3(Xsurvey,Ysurvey,Zsurvey,'*k')
h1 = plot3([Xsurvey(1) Xsurvey(1)], [Ysurvey(1) Ysurvey(1)], [zTest(1) Zsurvey(1)],'o-g');
h1 = plot3([Xsurvey(1)], [Ysurvey(1)], [Zsurvey(1)],'oc');
h1 = plot3([Xsurvey(end) Xsurvey(end)], [Ysurvey(end) Ysurvey(end)], [zTest(end) Zsurvey(end)],'o-g');
h2 = plot3(xiRunup,yiRunup,ziRunup,'*r');
xlim([min(cams.x,min(Xsurvey)), max(Xsurvey)+10])
ylim([ymin, ymax])
zlim([-1, 5])
xlabel('Cross-shore (m)')
ylabel('Alongshore (m)')
zlabel('Elevation (m)')
view([10 7.5])
title(['Madeira Beach - Line Number ',line_number,' - Survey Date: ',survey_date])
print('-dpng',['\\gs\stpetersburgfl-g\NACCH\Imagery\madbeach\runup\runupProjected\runupProjected_',camera_number,'_line',line_number,'.png'])
hgsave(       ['\\gs\stpetersburgfl-g\NACCH\Imagery\madbeach\runup\runupProjected\runupProjected_',camera_number,'_line',line_number,'.fig'])

figure
subplot(211)
plot(ziRunup,'*r')
hold on; grid on
plot(zTest,'*g')
xlim([0 length(ziRunup)])
ylabel('Runup line elevation (m)')
legend('Runup projected on beach', 'Runup verically interpolated','location','best')
title(['Madeira Beach - Line Number ',line_number,' - Survey Date: ',survey_date])

subplot(212)
plot(ziRunup-zTest)
grid on
xlim([0 length(ziRunup)])
xlabel('Pixel index of runup line')
ylabel('projection - interpolated (m)')
print('-dpng',['\\gs\stpetersburgfl-g\NACCH\Imagery\madbeach\runup\runupProjected\runupProfileProjectedVSInterp_',camera_number,'_line',line_number,'.png'])
hgsave(       ['\\gs\stpetersburgfl-g\NACCH\Imagery\madbeach\runup\runupProjected\runupProfileProjectedVSInterp_',camera_number,'_line',line_number,'.fig'])

%% now convert xiRunup, yiRunup to UV coords
% need to convert back to XYZ first...
clear XYZ
theta = ((90-47.753111137952047+90)); %angle of rotation, degrees
E0 = 3.230550667500000e+05; %origin, E
N0 = 3.075922214137756e+06; %origin, N

localOrigin = [E0,N0];
localAngle = 360-theta; %
directionFlag = 0;
Xin = xiRunup; Yin = yiRunup;
tmp=isnan(Xin);
Xin(tmp)=[];Yin(tmp)=[];ziRunup(tmp)=[];
[Xout Yout]= localTransformPoints(localOrigin,localAngle,directionFlag,Xin,Yin);

XYZ(:,1)=Xout;
XYZ(:,2)=Yout;
XYZ(:,3)=ziRunup;
tmp=isnan(XYZ(:,1));
XYZ(tmp,:)=[];

XYZworld = XYZ;
XYZlocal = [Xin' Yin' ziRunup'];

load(['\\gs\stpetersburgfl-g\NACCH\Imagery\madbeach\calibration\intrinsic\I2RGUS\20220125\',camera_number,'_madbeach_IOEOBest.mat'])

[UVd,flag] = xyz2DistUV(intrinsics,extrinsics,XYZ);
UVd = reshape(UVd,[],2);
s=size(XYZ(:,1));
Ud=(reshape(UVd(:,1),s(1),s(2)));
Vd=(reshape(UVd(:,2),s(1),s(2)));

% Round UVd coordinates so it cooresponds to matrix indicies in image I
Ud=round(Ud);
Vd=round(Vd);

% Utilize Flag to remove invalid points (those outside image, etc)
Ud(find(flag==0))=nan;
Vd(find(flag==0))=nan;

UVd = [Ud Vd];
tmp=isnan(UVd(:,1));
UVd(tmp,:)=[];
XYZworld(tmp,:) = [];
XYZlocal(tmp,:) = [];

save(['.\runup\runupProjected\',camera_number,'_runup_surveyline',line_number,'.mat'],'UVd','XYZworld','XYZlocal')
dlmwrite([camera_number,'_runup_surveyline',line_number,'.pix'],UVd,' ')
