% Transform UTM to a local "madbeach" coordinate system
%
% Input:
%   E, N = UTM coordinates
%
% Output:
%   X, Y = local coordinates
%
% Usage: [X,Y]=coordSys_madbeach(E,N)
%
% Written By:
%   Jenna Brown, 2017, USGS
%       based on Dave Thompson's madbeach_transformUTMToLocal.m
%--------------------------------------------------------------------------
function [X,Y]=coordSys_madbeach(E,N)

%from calibration extrinsics and intrinsics - 1/28/2022
theta = 132.246888862048; %angle of rotation, degrees
E0 = 323055.06675; %origin, E
N0 = 3075922.21413775; %origin, N

%[lon0, lat0] = UTM2ll(E0, N0, 17, 23);
%lon0 = -82.796102542950635;
%lat0 = 27.796216949206798;

[X,Y] = xyRotate(E,N,theta,E0,N0);

% %%%%%%testing rotatePoints() instead of xyRotate
% points = [E, N];
% origin = [E0, N0];
% XY = rotatePoints(points, theta, origin, 'degs');
% X = XY(:,1);
% Y = XY(:,2);


