function params = readparams(inputfile)

fid = fopen(inputfile,'r');
C = textscan(fid,'%s %s\n');
fclose(fid);

params = struct(C{1,1}{1},C{1,2}{1});

for i = 2:size(C{1,1},1)
    params.(C{1,1}{i}) = C{1,2}{i};
end

if ~endsWith(params.WS_folder,'/')
    params.WS_folder = [params.WS_folder '/'];
end

if ~endsWith(params.WS_data,'/')
    params.WS_data = [params.WS_data '/'];
end

if ~endsWith(params.WS_result,'/')
    params.WS_result = [params.WS_result '/'];
end

if ~endsWith(params.WS_browse,'/')
    params.WS_browse = [params.WS_browse '/'];
end

params.DTM_width = str2double(params.DTM_width);

params.LG_no_frequency = str2double(params.LG_no_frequency);
params.LG_no_orientation = str2double(params.LG_no_orientation);
params.LG_multiple_factor = str2double(params.LG_multiple_factor);

scaleStr = strsplit(params.CWT_scales,',');
params.CWT_scales = zeros(1,length(scaleStr));
for i = 1:length(scaleStr)
    params.CWT_scales(i) = str2num(scaleStr{i});
end
params.CWT_wavelet = lower(params.CWT_wavelet);

params.SF_dist = str2num(params.SF_dist);

params.PT_mintrace = str2num(params.PT_mintrace);
params.PT_minpts = str2num(params.PT_minpts);
params.PT_sigma = str2double(params.PT_sigma);
end