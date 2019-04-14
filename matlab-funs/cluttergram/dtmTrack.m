%% Open MOLA DTM, and interpolate it by using x,y coordinates
% Si-Ting Xiong
% Compiled on 11-03-2017


function z = dtmTrack(dtmfile,x,y)
[dtm,dtminfo] = readDTM(dtmfile);

%[x,y] = subTrackGeo(filepath, trackID, area);

z = dtmInterpxyz(dtm,ulx,uly,resx,resy,x,y);

end
