function procSHARAD(inputfile)
%% adding matlab packages
addpath('./matlab-funs/sharad_io/');
addpath('./matlab-funs/sharad_layers/');
addpath('./matlab-funs/utilities/');
addpath('./matlab-funs/cluttergram/');


%% define input files and folders
params = readparams(inputfile);

DT = datetime('now');
DateString = datestr(DT);
logfile = [params.WS_folder '/procSHARAD_' DateString '.log'];
fid = fopen(logfile,'w');

fprintf(fid,'ROI:%s\n',params.ROI_shapefile);
fprintf(fid,'Radargram list %s\n',params.WS_filelist);
fprintf(fid,'Data folder %s\n',params.WS_data);
fprintf(fid,'DTM file is %s\n',params.DTM);
fprintf(fid,'Results are outputed into %s\n',params.WS_result);

%% read geometry from shapfile
fprintf(fid,'Processing: reading ROI ...\n');
S = shaperead(params.ROI_shapefile);
if S.Geometry ~= 'Polygon'
    return
end
node = [S.X' S.Y'];
n      = size(node,1);
node = node(1:n-1,:);
n      = size(node,1);
cnect  = [(1:n-1)' (2:n)'; n 1];
fprintf(fid,'End\n');

%% read DTM...
[dtm,dtminfo] = readDTM(params.DTM);

%% loop for every radargram in the datapool
fileID = fopen(params.WS_filelist);
SHARADlist = textscan(fileID,'%s');
fclose(fileID);
SHARADlist = SHARADlist{1};

fprintf(fid,'Looping starting...\n');
if iscellstr(SHARADlist)
    for i = 507:length(SHARADlist)
        orbitid = SHARADlist{i};
        fprintf(fid,'####################################################\n');
        fprintf(fid,'Processing #%d: %s\n',i,orbitid);
        disp(orbitid);
        procOneOrbit(fid,params,orbitid,node,cnect,dtm, dtminfo);

        fprintf(fid,'Processing End #%d\n',i);
        fprintf(fid,'####################################################\n');
    end
end
fclose(fid);

end