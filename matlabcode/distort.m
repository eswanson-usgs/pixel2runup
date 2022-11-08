function [Ud, Vd] = distort(Uu, Vu, cam, ip)
%
% [Ud, Vd] = distort(Uu, Vu, cam, ip)
%
%  Routine to re-distort pixel coordinates for display on an image.
%  This version assumes that the image processor and camera data are stored
%  within structures (and can be retrieved by DBgetImageData or equiv).
%  Uu and Vu are columns of undistorted pixel coordinates, 
%  while Ud and Vd are their distorted equivalents.
%  Shorter calls are
% 	[Ud, Vd] = distort(UVu, cam, ip);
% 	UVd = distort(Uu, Vu, cam, ip);
% 	UVd = distort(UVu, cam, ip);
%
%  where UV is a matrix of 2 columns [U, V].
%
%  This version corrects a previous problem wherein image coordinate which 
%  were very far off-screen could be sometimes re-distorted back into the 
%  image.  This is done by omitting stopping the correction for points
%  whose distance from image center is beyond the turnaround in the UV
%  vs newUV function.
%

%  From several previous versions, by Holman, 7/16/98

% input error checking

global whichDistort;

if nargin == 3      % UV passed as matrix
    ip = cam;
    cam = Vu;
    if size(Uu,2) ~= 2
        Uu = Uu'
    end
    if size(Uu) ~= 2
        error('invalid arguments for distort ([Uu Vu] or Uu, Vu expected')
    end
    Vu = Uu(:,2);     % sort into columns
    Uu = Uu(:,1);
end
U = Uu(:);
V = Vu(:);
if length(U) ~= length(V)
    error('U and V must be the same length (distort)')
end

% here we have sorted everything out, cam, ip, UV, etc.
% new db: ip is basically useless and empty. cam has all the data, in a new format.
if (isfield(cam,'Drad') ...
	&& ~isempty(cam.K) ...
	&& (numel(cam.K)>1) ...
	&& (numel(cam.Drad)>2) )
	whichDistort = '2';
	[Ud, Vd] = distort2( U, V, cam, ip );
	return;
end;

whichDistort = '1';

%  Compute the maximum d for which you will allow distortion.  Beyond this
%  value distortion can falsely bring points back onto the screen.
%  Calculation is just based on finding max in newU vs U analytically

if cam.D1 ~= 0
	d2max = -(1 + cam.D2) / (3 * cam.D1);
else
	d2max = 1e8;	% some large number
end

%  Now carry out the distortion calculation.  Default is to return original points

x = (U / ip.lx) - ip.U0;
y = (V / ip.ly) - ip.V0;
d2 = x.^2 + y.^2;

good = d2 < d2max;
Ud = U;
Vd = V;

scale = 1 + cam.D2 + cam.D1*d2(good);
Ud(good) = ((x(good) .* scale) + ip.U0) * ip.lx;
Vd(good) = ((y(good) .* scale) + ip.V0) * ip.ly;

if nargout == 1
  Ud = [Ud Vd];
end

%
% Copyright by Oregon State University, 2002
% Developed through collaborative effort of the Argus Users Group
% For official use by the Argus Users Group or other licensed activities.
%
% $Id: distort.m 6 2016-02-11 00:46:00Z  $
%
% $Log: distort.m,v $
% Revision 1.7  2016/02/10 22:13:24  stanley
% better checks on K and Drad for validity
%
% Revision 1.6  2011/04/15 21:04:23  stanley
% verify non-empty K for new distort
%
% Revision 1.5  2008/01/22 23:19:37  stanley
% added test for newer db format, call distort2
%
% Revision 1.4  2006/08/11 22:36:37  stanley
% fixed camera/cam in 29 and 30
%
% Revision 1.3  2006/08/11 00:08:55  stanley
% fixed d2max calc by sign, from L.Clarke
%
% Revision 1.2  2006/04/04 18:46:42  stanley
% dmoraru -- remove unneeded sqrt
%
% Revision 1.1  2004/08/18 20:55:21  stanley
% Initial revision
%
%
%key geomtool internal 
%comment  
%