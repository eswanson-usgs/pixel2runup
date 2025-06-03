function [Ud, Vd] = distort2(Uu, Vu, cam, ip)
%
% version 2, for new database with cam containing IP
%
% [Ud, Vd] = distort(Uu, Vu, cam, ip)
%
%  Routine to re-distort pixel coordinates for display on an image.
%  This version assumes that the camera data are stored
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

if nargin == 2      % UV passed as matrix
    cam = Vu;
    if size(Uu,2) ~= 2
        Uu = Uu'
    end
    if size(Uu) ~= 2
        error('invalid arguments for undistort ([Uu Vu] or Uu, Vu expected')
    end
    Vu = Uu(:,2);     % sort into columns
    Uu = Uu(:,1);
end
U = Uu(:);
V = Vu(:);
if length(U) ~= length(V)
    error('U and V must be the same length (distort)')
end

%d2max = 1e8;	% some large number
% let's be smarter about our maximum radius that we allow
%  we have D_U0, D_V0, center of image. Maximum radius we could have is
%  close to sqrt(D_U0*D_V0) (not IS, because image center can be
%  left,below of real image "center". So, add 20%. Nah, 50%.
% oops, forgot that we are in K units, so must also include K

d2max = ((cam.D_U0/cam.K(1,1))^2 + (cam.D_V0/cam.K(2,2))^2) * 1.5 ;

%  Now carry out the distortion calculation.  Default is to return original points

x = (U - cam.D_U0) / cam.K(1,1);
y = (V - cam.D_V0) / cam.K(2,2);
%  note to John: K(2,2) is negative, but negatives cancel in line above and
%  calculation of Vd(good) below, so is ok. 
d2 = x.^2 + y.^2;
r = sqrt(d2);

good = d2 < d2max;
Ud = U;
Vd = V;

warning off MATLAB:divideByZero
scale = 1 + ( polyval( cam.Drad, r(good)) ./ r(good) );
warning on MATLAB:divideByZero

% is possible to have a radius of zero! at r=0, scale = 1 by definition.
scale(find(r(good)==0)) = 1;
Ud(good) = ((x(good) .* scale) * cam.K(1,1)) + cam.D_U0;
Vd(good) = ((y(good) .* scale) * cam.K(2,2)) + cam.D_V0;

% remove off-screen answers -- cannot remove, some things are counting
% on getting the same number of answers as questions. Sigh.
%Ud(find(Ud<1)) = 1;
%Vd(find(Vd<1)) = 1;
%Ud(find(Ud>ip.width)) = ip.width;
%Vd(find(Vd>ip.height)) = ip.height;

if nargout == 1
  Ud = [Ud Vd];
end

%
% Copyright by Oregon State University, 2002
% Developed through collaborative effort of the Argus Users Group
% For official use by the Argus Users Group or other licensed activities.
%
% $Id: distort2.m 6 2016-02-11 00:46:00Z  $
%
% $Log: distort2.m,v $
% Revision 1.3  2011/04/15 21:03:11  stanley
% remove clipping so offscreen will be off screen
%
% Revision 1.2  2010/11/17 19:03:02  stanley
% added tests for zero radius, UV off image.
%
% Revision 1.1  2008/01/22 23:21:43  stanley
% Initial revision
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