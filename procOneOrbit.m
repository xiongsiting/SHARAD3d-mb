function [flag,in_ind] = procOneOrbit(fid,params,orbitid, node, cnect, dtm, dtminfo)
                                 
LIGHTSPEED = 299792458; TSPACE = 0.0375e-6;

prefix = ['s_' orbitid];
fprecision = '%14.4f';

% Open GEO TAB file
fprintf(fid,'Subprocessing: reading geoinfo from .tab file for %s\n',orbitid);
pjgeofile = [prefix '_geom-proj.tab'];
pjgeofile = [params.WS_data 'proj-geom/' pjgeofile];
fileID = fopen(pjgeofile);
G = textscan(fileID,'%f %q %f %f %f %f %f %f %f %f %f %f','Delimiter',',');
                    % col, time, lat,lon,mrad,srad,rvel,tvel,sza,phase
fclose(fileID);
p = [G{11} G{12}];
fprintf(fid,'Subprocessing End.\n');

% Call inpoly
fprintf(fid,'Subprocessing: extracing indices within study area...\n');
% find which footprints are within the ROI
in = inpoly(p,node,cnect);
in_ind = find(in==1);
if isempty(in_ind) || length(in_ind) < params.PT_mintrace
    fprintf(fid,'Warning: %s has %d footprint locating inside study area\n',orbitid,length(in_ind));
    fprintf(fid,'Subprocessing End.\n');
    return;
end
fprintf(fid,'Subprocessing: %s has %d footprint locating inside study area\n',orbitid,length(in_ind));
fprintf(fid,'TrackID column No.: %s %d %d\n',orbitid, in_ind(1),in_ind(end));
fprintf(fid,'Subprocessing End.\n');

% Convert lat & lon to xy coordinates
geoinfo = [G{3} G{4} G{11} G{12} G{5}*1000 G{6}*1000];
geoinfo = geoinfo(in_ind,:);
z = dtmInterpxyz(dtm,dtminfo,geoinfo(:,3),geoinfo(:,4));

% Output the ll & xy to a ASCII file
locfile = [params.WS_result prefix '_llxyhh.txt'];
fprintf(fid,'Outputing: outputing the coordiates (lat,lon,x,y,Ha,Hsat) to %s\n',locfile);
dlmwrite(locfile, geoinfo,'precision',fprecision);
heifile = [params.WS_result prefix '_dtmh.txt'];

fprintf(fid,'Outputing: outputing DTM height to %s\n',heifile);
dlmwrite(heifile,z,'precision',fprecision);
% Subtrack the geofile and imgfile and output to a location

% Read IMG file and subtrack image
fprintf(fid,'Subprocessing: reading radargram from .img file for %s\n',orbitid);
im = readImg(params.WS_data,orbitid); %% has this changed to dB or not
if size(im,1) ~= 3600
    return;
end
inim = im(:,in_ind);
if ~isempty(find(inim == 0,1))
    inim(inim == 0) = 10e-12;
    dbim = 10*log10(inim);
else
    dbim = 10*log10(inim);
end
imgray = mat2gray(dbim);
fprintf(fid,'Subprocessing End.\n');

%%% check if this works for situations when the data is blank in the 'in' area
if isempty(find(imgray > 0, 1))  
    flag = 0;
    in_ind = [];
    return;
end


%%%%%%%%%%%%%%%%%%%%%%% Long time processing, be careful! %%%%%%%%%%%%%%%%%%%

% % % % %% Pick up surface echo from the radargram
% % % % surfecho = pickFirstReturn(D, 'mouginot');
% % % % cols = 1:length(in_ind);
% % % % indices = sub2ind(size(D),double(round(surfecho)),cols');
% % % % echoes = dbim(indices); %% using imgray or using original data???
% % % % % output the picked surface echo to outfolder
% % % % echofile = [outfolder prefix];
% % % % echofile = [echofile '_surfecho.txt'];
% % % % fprintf(fid,'Outputing: outputing the picked surface echoes to %s\n',echofile);
% % % % dlmwrite(echofile,[surfecho(rmline+1:size(surfecho,1)-rmline,:) echoes(rmline+1:size(echoes,1)-rmline,:)],'precision','%14.4f'); t1 = toc;
% % % % fprintf(fid,'Subprocessing End at %.6f.\n',t1);


%% Generating cluttergram using MOLA or HRSC DTMs
fprintf(fid,'Subprocessing: generating cluttergram from MOLA DTM...\n'); tic;
% resample the DTM according to the SHARAD track footprints
[dtmpatch,X,Y] = dtmPatch(dtm,dtminfo,geoinfo,params.DTM_width);
% Write out the dtm patch into a tif image
tiffile = [params.WS_result prefix '_dtmpatch.tif'];
t = Tiff(tiffile,'w');
dtmpatch = single(dtmpatch);
tagstruct.ImageLength = size(dtmpatch,1);
tagstruct.ImageWidth = size(dtmpatch,2);
tagstruct.Compression = Tiff.Compression.None;
tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample = 32;
tagstruct.SamplesPerPixel = 3;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software = 'MATLAB';
t.setTag(tagstruct);
t.write(single(cat(3,dtmpatch, X, Y)));
t.close();

% Generate the cluttergram from dtmpatch
[Ps,surfind,surfhei] = cluttergram(dtmpatch,X,Y,geoinfo, '1'); 
% Write out the cluttergram into binary or tif image
tiffile = [params.WS_result prefix '_clutim.tif'];
t = Tiff(tiffile,'w');
Ps = single(Ps);
tagstruct.ImageLength = size(Ps,1);
tagstruct.ImageWidth = size(Ps,2);
tagstruct.Compression = Tiff.Compression.None;
tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample = 32;
tagstruct.SamplesPerPixel = 1;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software = 'MATLAB';
t.setTag(tagstruct);
t.write(single(Ps*10e22));
t.close();
% output the surface index and height file 
clutfile = [params.WS_result prefix '_clutsrf.txt'];
dlmwrite(clutfile,[surfind' surfhei'],'precision',fprecision);
fprintf(fid,'Outputing: outputing the picked surface echoes to %s\n',clutfile);

t1 = toc;
fprintf(fid,'Subprocessing End at %.6f.\n',t1);

%% log-Gabor Filtering (optional)
fprintf(fid,'Subprocessing: Filtering radargram by using log-Gabor filters.\n');tic;
k = 2;
nscale = params.LG_no_frequency;
mult = params.LG_multiple_factor;
norient = params.LG_no_orientation;
softness = 1;
D = imgray;
imgray = noisecomp(D,k,nscale,mult,norient,softness); t1 = toc;
fprintf(fid,'Subprocessing End at %.6f.\n',t1);

%% Pick up surface echo from the radargram
fprintf(fid,'Subprocessing: picking out surface echoes...\n'); tic;
surfecho = size(1,size(D,2));
for i = 1:size(D,2)
    segment = D(surfind(i) - params.SF_dist:surfind(i) + params.SF_dist,i);
    [~, maxind] = max(segment);
    surfecho(i) = surfind(i) - params.SF_dist - 1 + maxind;
end
surfindex = sub2ind(size(D),surfecho,1:size(D,2));
surfinten = dbim(surfindex);
surfechofile = [params.WS_result prefix '_surfecho.txt'];
dlmwrite(surfechofile,[surfecho' surfinten'],'precision',fprecision);
fprintf(fid,'Outputing: outputing the picked surface echoes to %s\n',surfechofile);

%% Extract subsurface echoes by using the tracing layer method
fprintf(fid,'Subprocessing: picking out subsurface echoes...\n'); tic;
scaleRange = params.CWT_scales;
wavelet = params.CWT_wavelet;
tic; cwtfeatim = peakimcwt(imgray,scaleRange,wavelet,surfecho);toc;

figure; imagesc(imgray);hold on;
BW = cwtfeatim > 0;
BW2 = bwareaopen(BW, params.PT_minpts);
cwtfeatim = cwtfeatim.*BW2;
%plot(surfecho,'w-');
plot(surfind,'w-');
[x,y] = find(cwtfeatim > 0);
plot(y,x,'r.'); 

%% remove the clutter reflection
cwtfeatim = removeclutter(double(Ps),surfind, cwtfeatim, surfecho);


%% fillingin codes to write out the feature images
fprintf(fid,'Outputing: outputing the picked subsurface echoes to \n');
% BW = cwtfeatim > 0;
% BW2 = bwareaopen(BW, params.PT_minpts);
% cwtfeatim = cwtfeatim.*BW2;
[x,y] = find(cwtfeatim > 0);
if ~isempty(x)
plot(y,x,'g.'); 
ylim([min(surfind)-200 max(surfind) + 300]);
%axis equal;
outpngfile = [params.WS_browse orbitid];
saveas(gcf,outpngfile,'png');
close(gcf);close all;

pt = zeros(length(x),10);
for n = 1:size(pt,1)
    pt(n,1:4) = geoinfo(y(n),1:4);
    pt(n,5) = geoinfo(y(n),5) + (1800-x(n)) * TSPACE * LIGHTSPEED/2;
    pt(n,6) = (x(n) - surfind(y(n))) * TSPACE * LIGHTSPEED /2 /sqrt(params.PT_sigma);
	pt(n,10) = x(n) - surfind(y(n));
	pt(n,11) = dbim(surfind(y(n)),n);
	end
indices = sub2ind(size(cwtfeatim),x,y);
pt(:,7) = dbim(indices);
pt(:,8) = x;
pt(:,9) = y;

ptfile = [params.WS_result prefix '_subpt.txt'];
dlmwrite(ptfile,pt,'precision',fprecision);
fprintf(fid,'Outputing: outputing the picked subsurface echoes to %s\n',ptfile);
t1 = toc;
else
    fprintf(fid,'Warning: there is no point extracted!');
end
fprintf(fid,'Subprocessing End at %.6f.\n',t1);
end
