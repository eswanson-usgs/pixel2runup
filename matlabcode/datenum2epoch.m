function epoch = datenum2epoch(dnum)
% datenum2epoch.m  A function to convert Matlab date numbers to Epoch
%                  (Unix) time stamps.
%
% USAGE:        epoch = datenum2epoch(dnum);
%
%       WHERE:  dnum is a 1D array of Matlab date numbers (days since
%                    January 1, 0000 at midnight).
%               epoch is a 1D array of Epoch (Unix) time stamps (seconds
%                     since January 1, 1970 at midnight).
%
% HISTORY:
% C. Sullivan,  09/10/07,   version 1.0

epochTimeBase = datenum(1970,1,1,0,0,0);
epoch = (dnum - epochTimeBase)*3600*24;