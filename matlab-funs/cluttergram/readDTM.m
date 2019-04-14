function [dtm,dtminfo] = readDTM(filepath)

if nargin < 1 || isempty(filepath)
    filepath = './Data/megr_s_512_mola_proj.tif';
end

[dtm,R] = geotiffread(filepath);
[nx,ny] = size(dtm);

%% extend to be read from tif
ulx = R.XWorldLimits(1);
uly = R.YWorldLimits(2);
lrx = R.XWorldLimits(2);
lry = R.YWorldLimits(1);

resx = (lrx - ulx) / ny;
resy = (lry - uly) / nx;

dtm(dtm > 1e7) = nan;
dtm = double(dtm);

% %%%%% not used in this code %%%%%%%%%
% ullat = -73.1519;       % See from megr_s_512.lbl for first line lat & lon
% ullon = 315.0000;       % See from megr_s_512.lbl for first line lat & lon
% % This has an error when used to interpolate with SHARAD location of LAT &
% % LON
% % The fifth parameter should be -89.99999 to avoid Inf number
% [ulx,uly] = polarstereo_fwd(ullat,ullon,3396000,0,-89.999999,0);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% res = 0.115e3;      % Resolution of MOLA DTM over polar region
                    %See MAP_SCALE in megr_s_512.lbl
% resx = res; resy = -res;
% ulx = - resx * nx/2;
% uly = - resy * ny/2;

dtminfo = [ulx,uly,resx,resy];
end
