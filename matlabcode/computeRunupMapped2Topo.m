%Convert horizontal runup (x,y) to vertical runup (z) using lidar data
%
% This function:
% - transforms lidar postion data into horiz/vert angles relative to an origin
% - transforms runup postion data into horiz/vert angles relative to an origin
%       NOTE: origin is defined as mean camera position
%
%       lidar:    (x,y,z) ---> (dx,dy,dz) ---> (dr,phi,alpha)
%       runup:    (x,y,z) ---> (dx,dy,dz) ---> (dr,phi,alpha)
%
% - computes topo surfaces (Fx,Fy,Fz) that pass through lidar points (phi,alpha) 
%       with values of lidar data (x,y,z), such that:
%
%                   xTopo = Fx(phiTopo,alphaTopo)
%                   yTopo = Fy(phiTopo,alphaTopo)
%                   zTopo = Fz(phiTopo,alphaTopo)
%
% - evaluates topo interpolant surfaces (Fx, Fy, Fz) at runup points (phi,alpha)
%        to determine interpolated runup values
%
%                   xiRunup = Fx(phiRunup,alphaRunup)
%                   yiRunup = Fy(phiRunup,alphaRunup)
%                   ziRunup = Fz(phiRunup,alphaRunup)
%
% Inputs:
%   xRunup, yRunup, zRunup = runup data
%   xTopo, yTopo, zTopo = lidar data
%   xCamera, yCamera, zCamera = camera position (x,y,z)
%
% Outputs:
%   xiRunup, yiRunup, ziRunup = interpolated runup data (to lidar data)
%
% Usage: [xiRunup,yiRunup,ziRunup]=computeRunupMapped2Topo(xRunup,yRunup,zRunup,xTopo,yTopo,zTopo,xCamera,yCamera,zCamera);
%--------------------------------------------------------------------------
function [xiRunup,yiRunup,ziRunup]=computeRunupMapped2Topo(xRunup,yRunup,zRunup,xTopo,yTopo,zTopo,xCamera,yCamera,zCamera)

%*** Define Origin for Mapping ***
%*********************************

% !!!!! camera positions [t, x, y z] as of April 2015:
% c1=[1424400000,	32.6, 585.64, 43.81];
% c2=[1424448120,	32.56, 585.62, 43.1];
% c3=[1424400000,	34.4636, 583.6048, 42.8223];
% c4=[1424400000,	34.4395, 583.2903, 42.8271];
% c5=[1424448120,	34.67, 583.36, 43.1];
% c6=[1424400000,	32.721,	586.856, 43.152];

%c0 = [33.5757, 584.7285, 43.1352];                                          %mean camera position (of all cameras)
%c0 = [33.5158, 584.5388, 43.1399];                                          %mean camera positions for north-facing cameras (c1-c4)
%c0 = [0, 0, 18];


c0=[xCamera, yCamera, zCamera];

%*** Topo Data ***
%*********************************
%---remove bad values in topo data
good = find(~isnan(zTopo));                                                
xTopo = xTopo(good); 
yTopo = yTopo(good); 
zTopo = zTopo(good);
clear good

%---convert positions relative to origin
dx = xTopo-c0(1);
dy = yTopo-c0(2);
dz = zTopo-c0(3);

%---convert relative positions to angles
%%% Rob's code:
dr = sqrt(dx.*dx + dy.*dy);
alpha = atan(dz./dr);               %vertical angles (inclination)
phi = atan(dy./dx);                 %horizontal angles (azimuth)
%alpha = atan2(dz,dr);               %vertical angles (inclination)
%phi = atan2(dy,dx);                 %horizontal angles (azimuth)

%dr = sqrt(dx.*dx + dy.*dy + dz.*dz);
%alpha = acos(dz./dr);               %vertical angles (inclination)
%phi = atan(dy./dx);                 %horizontal angles (azimuth)

%---compute topo surfaces that pass through points (phi,alpha) with values of topo data
%       such that x = Fx(phi,alpha), y = Fy(phi,alpha), z = Fz(phi,alpha)
Fx = TriScatteredInterp(phi,alpha,xTopo);    %MATLAB: TriScatteredInterp will be removed in a future release.     
Fy = TriScatteredInterp(phi,alpha,yTopo);
Fz = TriScatteredInterp(phi,alpha,zTopo);
% Fx = scatteredInterpolant(phi,alpha,xTopo);       
% Fy = scatteredInterpolant(phi,alpha,yTopo);
% Fz = scatteredInterpolant(phi,alpha,zTopo);


%*** Runup Data ***
%*********************************
%---convert positions relative to origin
dxr = xRunup-c0(1);
dyr = yRunup-c0(2);
dzr = zRunup-c0(3);

%---convert relative positions to angles
%%% Rob's code:
drr = sqrt(dxr.*dxr + dyr.*dyr);
alphar = atan(dzr./drr);            %vertical angles (inclination)
phir = atan(dyr./dxr);              %horizontal angles (azimuth)
%alphar = atan2(dzr,drr);            %vertical angles (inclination)
%phir = atan2(dyr,dxr);              %horizontal angles (azimuth)

%drr = sqrt(dxr.*dxr + dyr.*dyr + dzr.*dzr);
%alphar = acos(dzr./drr);            %vertical angles (inclination)
%phir = atan(dyr./dxr);              %horizontal angles (azimuth)

%---evaluate topo interpolant surfaces (Fx, Fy, Fz) at runup points (phir,alphar)
%        to determine interpolated runup values
xiRunup = Fx(phir,alphar);                       
yiRunup = Fy(phir,alphar);
ziRunup = Fz(phir,alphar);

