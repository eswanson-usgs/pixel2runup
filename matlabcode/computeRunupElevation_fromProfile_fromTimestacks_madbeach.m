%Computes runup & TWL elevations from timestacks/timeseries for madbeach
%
%   - downloads tide and wave data for given day
%   - converts horizontal runup timeseries to runup elevation timeseries
%   - computes R2% statistics from runup timeseries
%   - computes TWL = tide + R2%
%   - computes R2% parameterization using Stockdon et al. (2006)
%
%   ==> saves Runup (timeseries, R2%, param, etc.) to "savePath" as "yyyymmdd_HHMM_runupyloc.mat"
%--------------------------------------------------------------------------
function computeRunupElevation_fromProfile_fromTimestacks_madbeach(runupFile,topoFile)

station='madbeach';

year = runupFile.date(1);
month = runupFile.date(2);
day = runupFile.date(3);

yearday=num2str(datenum(year,month,day)-datenum(year-1,12,31),'%3.3d');     %compute yearday (day of year, 0-365)
monthName=datestr([year,month,day,0,0,0],'mmm');                            %3-letter abreviation of the month
monthStr=num2str(month,'%2.2d');
dayStr=num2str(day,'%2.2d');

%--------------------------------------------------------------------------
% --- Tide & Wave Data ---

try
    [tides,waves]=getTidesWaves_madbeach([year,month,day],1,1);
catch
    [tides,waves]=getTidesWaves_madbeach([year,month,day],1,0);
end

%--------------------------------------------------------------------------
% --- Topo Data ---
survey.x=topoFile.profile.x(:);
survey.y=topoFile.profile.y(:);
survey.z=topoFile.profile.z(:);

%--------------------------------------------------------------------------
% --- Runup Data ---

runupData=dir([runupFile.path '\' yearday '_' monthName '.' dayStr '\*' station '.cx.r*.mat']); %person estimated runup
%---for each runup time series:
for ss=1:length(runupData)
    
    runupName = parseFilename(runupData(ss).name,'noLocal');

    yloc    = str2double(runupName.type(end-1:end)); %alongshore location of pixel transect
    runupHR = runupName.when(12:13);
    runupMN = runupName.when(15:16);
    
    % --- Runup Elevation ---
    %[Runup]=convertRunup_ArgusTimestack_madbeach(runupData(ss).name,survey);
    [Runup]=convertRunup_fromProfile_ArgusTimestack_madbeach(runupData(ss).name,survey);

    
    % --- R2% from Runup Statistics ---
    tide=interp1(tides.t,tides.measured,Runup.ts,'linear','extrap');
    Runup.tide      = mean(tide);
    
    %***de-tide runup timeseries
    Runup.eta      = Runup.zi - tide;
    
    %***compute setup
    Runup.setup      = mean(Runup.eta);
    
    %***compute statistics
    tRunup = matlab2Epoch(Runup.ts)-matlab2Epoch(Runup.ts(1));          %time series (seconds)
    
    good=find(isnan(Runup.eta)==0);
    Runup.eta=Runup.eta(good);
    tRunup=tRunup(good);
    
    % might run into issue w/ gap in tide data; use these lines to get
    % around it, but change #s to suit your gap in data
    %     tmp1=tides.t; tmp1(172:176)=[];
    %     tmp2=tides.measured; tmp1(172:176)=[];
    %     tides.measured(172:176)=interp1(tmp1,tmp2,tides.t(172:176))
    %     
    
    [~,~,Runup.R10,Runup.R5,Runup.R2,Runup.R13,~,~,~,~,~,~] = runupMaxValueStats(Runup.eta,tRunup);
    
    [~,~,Runup.S,~,~,~,Runup.Sinc,Runup.SIG,~,~] = runupFullStats(Runup.eta,tRunup,16);
    
    Runup.TWL = Runup.tide + Runup.R2;
    
    
    % --- X2% from Runup Statistics ---
    %***compute statistics
    
    [~,~,~,~,Runup.X2,~,~,~,~,~,~,~] = runupMaxValueStats(-Runup.x,tRunup);
    [~,~,~,~,Runup.Xi2,~,~,~,~,~,~,~] = runupMaxValueStats(-Runup.xi,tRunup);
    Runup.X2 = -Runup.X2;       %account for "runup" being in negative direction
    Runup.Xi2 = -Runup.Xi2;
   
    
    % --- R2% parameterization (Stockdon et al., 2006) ---
    
    %***local beach profile & slope
    %iloc=find(min(abs(survey.y(:,1)-(-Runup.yloc)))==abs(survey.y(:,1)-(-Runup.yloc)));
    %Profile(1,:)=topo.x(iloc,:);
    %Profile(2,:)=topo.z(iloc,:);
    
    Profile(1,:) = survey.x;
    Profile(2,:) = survey.z;
    
    [Runup.slope]      = calcRunupLocalBeachSlope([Runup.ts, Runup.zi],Profile);
    
    Runup.topo=topoFile;
    Runup.profile.x=Profile(1,:);
    Runup.profile.z=Profile(2,:);
    Runup.notes={'all elevation data (tides, bathy, runup) are NAVD88'};
    
    %***offshore wave conditions
    if ~isempty(waves)
        Hs=mean(interp1(waves.t,waves.H,Runup.ts,'linear','extrap'));
        Tp=mean(interp1(waves.t,waves.Tp,Runup.ts,'linear','extrap'));
        
        [Hs0,~,~] = reverse_shoaling(waves.depth,Hs,Tp); % uses THETA=0
        
        L0 = (9.81*Tp.^2)./(2.*pi);
        
        Runup.Hs  = Hs;
        Runup.Hs0 = Hs0;
        Runup.Tp  = Tp;
        Runup.I   = abs(Runup.slope)./sqrt(Hs0./L0);
        
        %***parameterization
        [Runup.param.R2, Runup.param.S, Runup.param.setup, Runup.param.Sinc, Runup.param.SIG, ~] = calcR2(Hs0,Tp,Runup.slope);
    else
        Runup.Hs  = NaN;
        Runup.Hs0 = NaN;
        Runup.Tp  = NaN;
        Runup.I   = NaN;
        Runup.param.R2    = NaN;
        Runup.param.S     = NaN;
        Runup.param.setup = NaN;
        Runup.param.Sinc  = NaN;
        Runup.param.SIG   = NaN;
    end
    
    %---Output:
    runupTemplate.yloc    = 0; %alongshore location of runup
    runupTemplate.t       = 0; %mean time as datenum
    runupTemplate.TWL     = 0;
    runupTemplate.R2      = 0;
    runupTemplate.R5      = 0;
    runupTemplate.R10     = 0;
    runupTemplate.R13     = 0;
    runupTemplate.setup   = 0;
    runupTemplate.S       = 0;
    runupTemplate.Sinc    = 0;
    runupTemplate.SIG     = 0;
    runupTemplate.slope   = 0;
    runupTemplate.tide    = 0;
    runupTemplate.Hs      = 0; %significant wave height
    runupTemplate.Hs0     = 0; %deep-water significant wave height
    runupTemplate.Tp      = 0; %peak wave period
    runupTemplate.I       = 0; %Irribarren number
    runupTemplate.param   = struct([]);
    runupTemplate.X2      = 0;
    runupTemplate.Xi2     = 0;
    
    %***timeseries output:
    runupTemplate.ts      = zeros(2048,1); %time as datenum
    runupTemplate.x       = zeros(2048,1); %cross-shore position
    runupTemplate.y       = zeros(2048,1); %alongshore position
    runupTemplate.xi      = zeros(2048,1); %interpolated cross-shore position
    runupTemplate.yi      = zeros(2048,1); %interpolated alongshore position
    runupTemplate.zi      = zeros(2048,1); %interpolated elevation
    runupTemplate.eta     = zeros(2048,1);
    
    runupTemplate.topo   = struct([]);
    runupTemplate.profile = struct([]);
    runupTemplate.notes   = '';
    
    Runup = orderfields(Runup, runupTemplate);
    
    %---Save:
    save([runupFile.path '\' yearday '_' monthName '.' dayStr '\' num2str(year) '' num2str(month,'%2.2d') '' num2str(day,'%2.2d') '_' runupHR '' runupMN '_runup' num2str(yloc) '_fromProfile.mat'],'Runup','-v7.3');

    fprintf(['Done Saving ',runupName.when,'\n'])
   
    clear runupName yloc runupHR runupMN bathyIndiv* Runup RunupIndiv tide tRunup iloc Profile* Hs Hs0 Tp L0 good
end