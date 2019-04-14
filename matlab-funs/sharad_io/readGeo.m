function [lat,lon,Ha,Hsat,time] = readGeo(filepath,trackID)
C = 299792458;
filename = ['s_',trackID];
filename = [filepath,filename];
%% Open GEO TAB file
fileID = fopen([filename,'_geom.tab']);
% col, time, lat,lon,mrad,srad,rvel,tvel,sza,phase
G = textscan(fileID,'%d %q %f %f %f %f %f %f %f %f','Delimiter',',');
fclose(fileID);
% Get lat and lon from GEO TAB file
lat = G{3};
lon = G{4};
time = 2*(G{6} - G{5})*1000/C;
Ha = G{5}*1000;
Hsat = G{6}*1000;
end