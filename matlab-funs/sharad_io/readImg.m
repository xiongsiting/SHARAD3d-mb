%% Read the radargram into Matlab, return radargram which in unit of dB
% Si-Ting Xiong
% Compiled on 11-03-2017

function D = readImg(filepath,trackID)

filename = ['s_',trackID];
filename = [filepath,filename];
%% Open IMG file with parameter read from lbl
% D = imread('s_00294501_rgram.tif');
% [nx,ny]=size(D);

fid = fopen([filename,'_rgram.lbl']);
line = fgetl(fid);
while ischar(line)
%     disp(line);
    line = line(~isspace(line));
    line = strsplit(line,'=');
    if size(line,2) == 2
        if strcmp(line{1},'LINES')
            ny = str2num(line{2});
        else
            if strcmp(line{1},'LINE_SAMPLES')
                nx = str2num(line{2});
            else
                if strcmp(line{1},'SAMPLE_BITS')
                    nbytes = str2num(line{2});
                end
            end
        end
    end
    line = fgetl(fid);
end
fclose(fid);
if nbytes == 32
    type = 'single';
else
    if nbytes == 16
        type = 'int';
    end
end


fid = fopen([filename,'_rgram.img']);
D = fread(fid,[nx,ny],type);
fclose(fid);
D = D';
%D = 10*log10(D);

end
