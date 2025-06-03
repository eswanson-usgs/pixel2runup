E=[]
N=[]
elevation=[]
for ii=1:length(files)
if strcmp(files(ii).name(7:9),'raw')
else
data=dlmread(files(ii).name,' ');
E=[E;data(:,1)];
N=[N;data(:,2)];
elevation=[elevation;data(:,3)];
end
end

scatter(E,N,20,elevation,'filled')
colorbar
box on
grid on

[X,Y]=coordSys_madbeach(E,N)

figure
scatter(X,Y,20,elevation,'filled')
grid on