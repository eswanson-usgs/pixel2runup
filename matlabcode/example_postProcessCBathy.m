clear,close all

dateRange = [2017,6,1;
             2017,6,30];

%--------------------------------------------------------------------------
days = [datenum(dateRange(1,1),dateRange(1,2),dateRange(1,3)):datenum(0,0,1):datenum(dateRange(2,1),dateRange(2,2),dateRange(2,3))];

cB0 = 0;

for dd = 1:length(days)
    
    D = datevec(days(dd));
    year  = D(1);
    month = D(2);
    day   = D(3);
    clear D
    
    yearday=datenum(year,month,day)-datenum(year-1,12,31);
    yearday=num2str(yearday,'%03d');
    monthName=datestr([year,month,day,0,0,0],'mmm');
    dayStr=num2str(day,'%2.2d');
    dirName = [yearday '_' monthName '.' dayStr];
    
    %---Kalman filtered cBathy  
    cBOutputPn = '\\gs\StPetersburgFL-G\NACCH\Imagery\madbeach\testing\cBathyUpdate_kalmanFilter';
    
    H = 1; %---assume wave height is 1 m if absence of better data
    
    
    cBInputPn = ['\\gs\StPetersburgFL-G\NACCH\Imagery\madbeach\testing\cBathyUpdate\' dirName];
    C = dir([DIR '\*cBathy*.mat']);
    
    if cB0==0
        bathy = 
        
    
    bathy = cBathyResults(1);
    for ii = 2:length(cBathyResults)
        priorBathy = bathy;
        bathy = cBathyResults(ii);
        bathy = KalmanFilterBathyNotCIL(priorBathy, bathy, H);
        eval(['save ' cBOutputPn,filesep,'test_',num2str(ii,'%2.2d'),'.mat', ' bathy'])
    end
    
    fns = dir([cBOutputPn,filesep,'test*.mat']);
    
    for i = 1: length(fns)
        load([cBOutputPn, filesep, fns(i).name])
        %figure;
        plotBathyCollectKalman(bathy)
        subplot(1,2,1);
        caxis([0 2]);
        pause(2)
    end
end