function [dtmpatch,X,Y] = dtmPatch(dtm,dtminfo,geoinfo,width)
%moladtmfile = './DTM/MOLA/megr_s_512.tif';
%[moladtm,moladtminfo] = readDTM(moladtmfile);

% Updated on 13-10-2017 adapted from previous version for package of SPLD-PL
resx = dtminfo(3);
resy = dtminfo(4);
ulx = dtminfo(1);
uly = dtminfo(2);

x = geoinfo(:,3);
y = geoinfo(:,4);

%% Calculate the four corners of the subarea on DTM used to simulate the cluttergram
kt = (y(end)-y(1))/(x(end)-x(1));
k = -1/kt;

xhalf = width * cos(abs(atan(k)));
yhalf = width * sin(abs(atan(k)));
xpixels = ceil(2 * xhalf / resx);
ypixels = ceil(2 * yhalf / -resy);

acrosspixels = max(xpixels,ypixels);

X = zeros(acrosspixels,length(x));
Y = zeros(acrosspixels,length(x));
dtmpatch = zeros(acrosspixels,length(x));
for i = 1:length(x)
    if acrosspixels == xpixels
        X(:,i) = x(i)-abs(xhalf):resx: x(i) + abs(xhalf);
        Y(:,i) = k * X(:,i) + y(i) - k * x(i);
    else
        Y(:,i) = y(i)+abs(yhalf):resy:y(i) - abs(yhalf);
        X(:,i) = (Y(:,i) - y(i) + k * x(i))/ k;
    end
    dtmpatch(:,i) = dtmInterpxyz(dtm,dtminfo,X(:,i),Y(:,i));
    %indices = find(dtmpatch(:,i) == 3363232,1);
    %if ~isempty(indices)
    %    dtmpatch(:,i) = dtmInterpxyz(moladtm,moladtminfo,X(:,i),Y(:,i));
    %end
end


