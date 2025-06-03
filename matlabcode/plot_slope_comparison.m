%% plot comparisons between beach slopes

% define beach slopes
meas_fs   = 0.05566; % min elevation to MHW of 0.198m
meas_bs_1 = 0.0378; %MHW(0.198) to "toe" @ 1.198m
meas_bs_2 = 0.0415; %MHW(0.198) to "toe" @ 1.54m
model     = 0.0376; % model slope for camera location (as of new grid ~April 9, 2017)


%pull in waves
year  = '2017';
month = '06';
day   = '1';
% for d = 1:31
%     day = num2str(d,'%02d')

% buoys
% offshore (St. Pete) = 42099
% Egmont (Key) = 42098
buoyname = 'offshore';%
buoynum  = '42099'; %

load(['\\gs\StPetersburgFL-G\NACCH\Imagery\madbeach\waves\waves_',year,month,day,'_NDBC',buoynum,'_',buoyname,'.mat'])

slopes = [meas_fs;meas_bs_1;meas_bs_2;model];
for ii = 1:length(slopes)
slope = ones(size(waves.H))*slopes(ii);

% now find R2
[R2(:,ii),S(:,ii),setup(:,ii), Sinc(:,ii), SIG(:,ii), ir(:,ii), R16(:,ii)] = calcR2(waves.H,waves.Tp,slope,0);

end


figure(1);clf
plot(waves.t,R2(:,1:3))
hold on
plot(waves.t,R2(:,4),'--')
datetick('x','mm/dd-HH','keepticks')
title(['Runup Compare Slopes (calcR2) - ',year,month,day])
xlabel('date-time')
ylabel('runup elevation')
legend(['meas fs    = ',num2str(meas_fs)],['meas bs1 = ',num2str(meas_bs_1)],['meas bs2 = ',num2str(meas_bs_2)],...
    ['model       = ',num2str(model)])

% saveas(gcf,['slope_compare_',year,month,day],'png')

% pause(1)
% end