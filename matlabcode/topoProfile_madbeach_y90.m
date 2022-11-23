%madbeach - load topo profile according to date and line number
%
%Input:
%   topoFile.date
%           .lineno
%           .path
%
%Output:
%   topoFile.profile
%--------------------------------------------------------------------------
function [topoFile] = topoProfile_madbeach_y90(topoFile)

%--- Hypack line60 is y=-90m runup timestack transect location
if length(topoFile.lineno)==1
    load([topoFile.path '\line' num2str(topoFile.lineno) '.xyz']);
    eval(['ENz = line' num2str(topoFile.lineno) ';'])

    %---convert to local coordinate system
[line.x,line.y]=coordSys_madbeach(ENz(:,1),ENz(:,2));
line.z = ENz(:,3);

[profile.x,inds]=sort(line.x);
profile.y=line.y(inds);
profile.z=line.z(inds);

else
    % line first line file
    load([topoFile.path '\line' num2str(topoFile.lineno(1)) '.xyz']);
    eval(['ENz1 = line' num2str(topoFile.lineno(1)) ';'])
    % line second line file
    load([topoFile.path '\line' num2str(topoFile.lineno(2)) '.xyz']);
    eval(['ENz2 = line' num2str(topoFile.lineno(2)) ';'])
    
    %---convert to local coordinate system
    [line1.x,line1.y]=coordSys_madbeach(ENz1(:,1),ENz1(:,2));
    line1.z = ENz1(:,3);
    [profile1.x,inds]=sort(line1.x);
    profile1.y=line1.y(inds);
    profile1.z=line1.z(inds);
    
    [line2.x,line2.y]=coordSys_madbeach(ENz2(:,1),ENz2(:,2));
    line2.z = ENz2(:,3);
    [profile2.x,inds]=sort(line2.x);
    profile2.y=line2.y(inds);
    profile2.z=line2.z(inds);
    
    %interpolate line48 to line49
    profile1.z_new = interp1(profile1.x,profile1.z,profile2.x);
    
    profile.x=profile2.x;
    profile.y=profile2.y;
    profile.z=mean([profile1.z_new profile2.z],2);
%     keyboard
    %***write something to load both lines and compute the mean
end

% figure;
% plot(profile.x,profile.z,'k','linewidth',2);
% set(gca,'xdir','reverse');
% xlabel('cross-shore (m)');
% ylabel('elevation (m)');

topoFile.profile = profile;