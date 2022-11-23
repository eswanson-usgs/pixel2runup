%% plot runup line on snap (or other) image
%
cd('\\gs\stpetersburgfl-g\NACCH\Imagery\madbeach\')
% load snapshot images
c1Snap = imread('\\gs\stpetersburgfl-g\Coastal_Change_Hazards\Archive\Data\2022\2022-305-FA\madbeach\1643898600.c1.snap.jpg');
c2Snap = imread('\\gs\stpetersburgfl-g\Coastal_Change_Hazards\Archive\Data\2022\2022-305-FA\madbeach\1643898600.c2.snap.jpg');

c1_020 = dlmread('\\gs\stpetersburgfl-g\NACCH\Imagery\madbeach\runup\madbeach-runup\Runup_UV_coords\UVc1_runup20.txt');
c1_090 = dlmread('\\gs\stpetersburgfl-g\NACCH\Imagery\madbeach\runup\madbeach-runup\Runup_UV_coords\UV_runup90.txt');
c1_150 = dlmread('\\gs\stpetersburgfl-g\NACCH\Imagery\madbeach\runup\madbeach-runup\Runup_UV_coords\UV_runup150.txt');
c1_250 = dlmread('\\gs\stpetersburgfl-g\NACCH\Imagery\madbeach\runup\madbeach-runup\Runup_UV_coords\UV_runup250.txt');

c1_pix_old1 = dlmread('\\gs\stpetersburgfl-g\NACCH\Imagery\madbeach\runup\pix\pixcamera1madbeach_20220729_90.pix');
c1_pix_old2 = dlmread('\\gs\stpetersburgfl-g\NACCH\Imagery\madbeach\runup\pix\pixcamera1madbeach_20220801_90_150.pix');
c2_pix_old1 = dlmread('\\gs\stpetersburgfl-g\NACCH\Imagery\madbeach\runup\pix\pixcamera2madbeach_20220729_393.pix');

figure(1);clf
i=imagesc(c1Snap);
hold on;box on
a=plot(c1_090(:,1),c1_090(:,2),'x-');
b=plot(c1_020(:,1),c1_020(:,2),'x-');
c=plot(c1_150(:,1),c1_150(:,2),'x-');
d=plot(c1_250(:,1),c1_250(:,2),'x-');

e=plot(c1_pix_old1(:,1),c1_pix_old1(:,2),'x-');
f=plot(c1_pix_old2(:,1),c1_pix_old2(:,2),'x-');
uistack(f,'bottom');uistack(e,'bottom');uistack(i,'bottom');
% set(gca,'ydir','reverse') % if only plotting lines

figure(2);clf
imagesc(c2Snap)
hold on;box on
e=plot(c2_pix_old1(:,1),c2_pix_old1(:,2),'x-');


