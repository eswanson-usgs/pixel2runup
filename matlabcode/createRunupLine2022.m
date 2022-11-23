%% create runup line
% use pixel toolbox
addpath('C:\Users\eswanson\OneDrive - DOI\Documents\CIRN\PIXel-Toolbox')
addpath('C:\Users\eswanson\OneDrive - DOI\Documents\CIRN\Support-Routines')

%% c1
PIXForget;
stationName = 'madbeachc1';
PIXSetStation(stationName)
zmsl = 0.0;
% create runup timestacks
% create new runup stacks at locations closer to cameras and have overlap
% between c1 and c2. perhaps create 2 new stack lines one just south of
% first groin (probably less interference by chairs/umbrellas) and one
% further south that has larger overlap between images (but might have lots
% of interference from chairs/umbrellas)
% '-20' and '-40' are just an initial guess, will have to experiment to
% find 'optimal' locations for new  lines
y_runup_locs=[-20 -40 -90 -150 -250];
xshore = [150]; xdune = [0];
rot =0; %rotation
dx = 0.1;
for i = 1:length(y_runup_locs)
[xr, yr, zr] = runupPIXArray(xshore,xdune,y_runup_locs(i),zmsl,rot);
instName = ['runup' num2str(fix(y_runup_locs(i)))];
idr(i) = PIXCreateInstrument(instName,'runup',PIXFixedZ+PIXUniqueOnly);
PIXAddPoints(idr(i), xr, yr, zr, instName );
end

% create a large matrix coverage for cBathy
x1 = 30; x2 = 300; dx = 3;
y1 = -320; y2 = -20; dy = 5;
name = 'cBathyArray';
% idb=PIXCreateInstrument('mBW','matrix', PIXFixedZ+PIXDeBayerStack);
idc = PIXCreateInstrument('mBW','matrix',PIXInterpUV);
PIXAddMatrix(idc,x1, y1, x2, y2, zmsl, dx, dy, name);

pid = PIXCreatePackage('madbeachc1', [idr idc] );

epoch = datenum2epoch(datenum(2022,01,28));
tide = 0;
r = PIXCreateR( pid, epoch, tide, 'fov');

% r.cams = input camera calibration parameters
% r.geoms = input camera geometry parameters
% r.ips = input image info

cams = {'c1', 'c2'};

r = PIXParameterizeR( r, cams, geoms, ips );
r = PIXRebuildCollect( r );

